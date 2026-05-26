import numpy as np
import argparse
import sys
import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, SASA, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

if not hasattr(np, 'int'):
    np.int = int

# 忽略非关键警告
warnings.simplefilter('ignore', BiopythonWarning)

class MultiChainImmuneScanner:
    def __init__(self, pdb_file, virus_chain_ids, receptor_chain_ids, ligand_name=None):
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('complex', pdb_file)
        self.model = self.structure[0]
        
        # 将输入字符串转换为列表 (e.g., "ABC" -> ['A', 'B', 'C'])
        self.virus_chain_ids = list(virus_chain_ids)
        self.receptor_chain_ids = list(receptor_chain_ids)
        self.ligand_name = ligand_name.upper() if ligand_name else None
        
        # 参数设定 parameter
        self.PROBE_SMALL = 1.4   # 水分子 water molecular
        self.PROBE_LARGE = 10.0  # 抗体探针 antibody probe
        self.FAB_RADIUS = 45.0   # 虚拟 Fab 半径 Fab radius
        
        # 数据容器
        self.virus_atoms = []
        self.receptor_atoms = []
        self.virus_residues = []
        
        self._parse_structure()

    def _parse_structure(self):
        #"""解析结构，提取病毒和受体原子"""
        print(f"Parsing Structure...")
        print(f"  - Virus Chains: {self.virus_chain_ids}")
        print(f"  - Receptor Chains: {self.receptor_chain_ids}")
        if self.ligand_name:
            print(f"  - Target Ligand: {self.ligand_name}")

        found_virus_chains = set()
        
        for chain in self.model:
            # --- 1. 处理病毒链 (多链) ---
            if chain.id in self.virus_chain_ids:
                found_virus_chains.add(chain.id)
                for residue in chain:
                    if is_aa(residue, standard=True):
                        self.virus_residues.append(residue)
                        for atom in residue:
                            self.virus_atoms.append(atom)
            
            # --- 2. 处理受体链 (配体逻辑) ---
            if self.ligand_name:
                 for residue in chain:
                    #print(residue.get_resname().strip().upper())
                    if residue.get_resname().strip().upper() == self.ligand_name:
                       for atom in residue:
                           self.receptor_atoms.append(atom)

            # --- 3. 处理受体链 (蛋白受体) ---
            if chain.id in self.receptor_chain_ids:
                for residue in chain:
                    if not self.ligand_name:
                        for atom in residue:
                            self.receptor_atoms.append(atom)

        # --- 错误检查 ---
        missing_virus = set(self.virus_chain_ids) - found_virus_chains
        if missing_virus:
            print(f"Warning: Virus chains {missing_virus} not found in PDB.")
        
        if not self.virus_atoms:
            print("Error: No virus atoms found.")
            sys.exit(1)
        
        if not self.receptor_atoms:
            print("Error: No receptor/ligand atoms found.")
            sys.exit(1)
            
        print(f"  - Loaded {len(self.virus_residues)} virus residues.")
        print(f"  - Loaded {len(self.receptor_atoms)} receptor atoms.")

    def calculate_sasa(self):
        #"""
        #计算 SASA。
        #关键策略：暂时将受体链从 Model 中 Detach (移除)，
        #保留完整的病毒复合物 (如三聚体) 进行计算，
        #这样能保证病毒内部界面 (Interface) 是掩埋的，而受体结合面是暴露的。
        #"""
        
        # 1. 暂时移除受体链 temporal remove of the receptor chains
        detached_receptor_chains = []
        for cid in self.receptor_chain_ids:
            if cid in self.model:
                chain_obj = self.model[cid]
                self.model.detach_child(cid)
                detached_receptor_chains.append(chain_obj)
        
        # 此时 self.model 只包含病毒链 (以及其他无关链)，
        # 病毒各链之间的相互遮挡会被正确计算。
        
        sr = SASA.ShrakeRupley()
        sr.search_points = 100
        
        # --- Round 1: Small Probe (Water) ---
        sr.compute(self.model, level='R')
        sasa_small_map = {res: res.sasa for res in self.virus_residues}
        
        # --- Round 2: Large Probe (Antibody) ---
        # 同样使用半径欺骗法
        old_radii = {}
        delta_radius = self.PROBE_LARGE - self.PROBE_SMALL
        
        for atom in self.virus_atoms:
            old_radii[atom] = atom.xtra.get('radius', 1.5)
            atom.xtra['radius'] = old_radii[atom] + delta_radius
            
        sr.compute(self.model, level='R')
        sasa_large_map = {res: res.sasa for res in self.virus_residues}
        
        # 恢复原子半径
        for atom in self.virus_atoms:
            atom.xtra['radius'] = old_radii[atom]
            # 恢复 residue.sasa 为小探针的值 (可选，为了数据一致性)
            res = atom.get_parent()
            if res in sasa_small_map:
                res.sasa = sasa_small_map[res]
        
        # 3. 将受体链装回去 (虽然碰撞计算用的是 coordinates list，但恢复结构是好习惯)
        for chain_obj in detached_receptor_chains:
            self.model.add(chain_obj)
            
        return sasa_small_map, sasa_large_map

    def calculate_clash(self):
        #"""计算虚拟抗体碰撞 (自适应受体尺寸)"""
        receptor_ns = NeighborSearch(self.receptor_atoms)
        clash_results = {}
        
        # 如果是小配体模式 (比如流感 SIA) -> 使用“抗体足迹覆盖”模型
        if self.ligand_name:
            # 假设抗体结合在你身上，抗体的身躯能覆盖周围 15A 的范围 (Fab的横向半径) 
            ANTIBODY_FOOTPRINT_RADIUS = self.FAB_RADIUS 
            
            for res in self.virus_residues:
                try:
                    res_center = res['CA'].get_coord()
                except KeyError:
                    atom_coords = [a.get_coord() for a in res]
                    res_center = np.mean(atom_coords, axis=0)
                
                # 直接搜索残基周围 15A 内有没有配体原子
                # 越近，说明抗体盖住配体的概率越高，Clash 分数越高
                nearby_ligand_atoms = receptor_ns.search(res_center, ANTIBODY_FOOTPRINT_RADIUS, level='A')
                
                if nearby_ligand_atoms:
                    # 计算到最近配体原子的距离
                    distances = [np.linalg.norm(res_center - atom.get_coord()) for atom in nearby_ligand_atoms]
                    min_dist = min(distances)
                    
                    # 距离越近，Clash 权重越高 (反比关系)
                    # 比如距离 3A 分数很高，距离 14A 分数很低
                    clash_score = max(0, int((ANTIBODY_FOOTPRINT_RADIUS - min_dist) * 10))
                else:
                    clash_score = 0
                    
                clash_results[res] = clash_score

        # 如果是大蛋白模式 (比如新冠 ACE2) -> 使用原来的“法向量外推大球”模型
        else:
            coords = [atom.get_coord() for atom in self.virus_atoms]
            virus_com = np.mean(coords, axis=0)
            
            for res in self.virus_residues:
                try:
                    res_center = res['CA'].get_coord()
                except KeyError:
                    atom_coords = [a.get_coord() for a in res]
                    res_center = np.mean(atom_coords, axis=0)
                
                normal_vec = res_center - virus_com
                norm = np.linalg.norm(normal_vec)
                if norm == 0: 
                    clash_results[res] = 0
                    continue
                unit_normal = normal_vec / norm
                
                sphere_center = res_center + (unit_normal * (self.FAB_RADIUS + 2.0))
                nearby_atoms = receptor_ns.search(sphere_center, self.FAB_RADIUS, level='A')
                clash_results[res] = len(nearby_atoms)
                
        return clash_results

    def run_analysis(self):
        print("\n--- Starting Analysis ---")
        sasa_small, sasa_large = self.calculate_sasa()
        clash_scores = self.calculate_clash()
        
        print("\nResID\tChain\tAA\tSASA_Sm\tSASA_Lg\tClash\tCategory")
        print("-" * 90)
        
        results = []
        
        for res in self.virus_residues:
            rid = res.get_id()[1]
            chain = res.get_parent().id
            rname = res.get_resname()
            
            ss = sasa_small.get(res, 0)
            sl = sasa_large.get(res, 0)
            clash = clash_scores.get(res, 0)
            
            # --- 自动分类逻辑 ---
            category = "Buried"
            
            # 只有当物理表面暴露时才讨论 (SASA_Small > 10)
            if ss > 10.0: 
                if sl < 5.0:
                    # 只有水能进，抗体进不去
                    category = "Inaccessible_Pocket"
                else:
                    # 抗体能碰到
                    if clash > 0:
                        category = "Competes_Receptor" # 黄金区域
                    else:
                        category = "Non_Neutralizing"
            
            # 简单的输出过滤：只打印有趣的位点 (可选)
            # if category != "Buried":
            print(f"{rid}\t{chain}\t{rname}\t{ss:.1f}\t{sl:.1f}\t{clash}\t{category}")
            
            results.append({
                'chain': chain, 'id': rid, 'aa': rname, 
                'sasa_small': ss, 'sasa_large': sl, 
                'clash': clash, 'category': category
            })
            
        return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multi-Chain Viral Immune Escape Scanner")
    parser.add_argument("--pdb", required=True, help="Path to PDB file")
    
    # 这里的 help 更新了：支持输入字符串，如 'ABC'
    parser.add_argument("--virus", required=True, help="Chain IDs of Virus (e.g., 'ABC' for trimer, 'A' for monomer)")
    parser.add_argument("--receptor", required=True, help="Chain IDs of Receptor (e.g., 'D', or 'HL' for antibody)")
    parser.add_argument("--ligand", required=False, help="[Optional] Ligand Name (e.g., SIA). If ignored, assumes Protein Receptor.")
    
    args = parser.parse_args()
    
    #SARS_CoV-2
    #pdb_file = "Q://KeyAntigenicSite/ViralStruct/SARS2/PDB/CoV/7DX5.pdb" 
    #Kvirus_chains='ABC'      # 病毒部分
    #Kreceptor_chains='D'    # 受体部分
    #scanner = MultiChainImmuneScanner(pdb_file, Kvirus_chains, Kreceptor_chains)

    #H3N2
    #pdb_file = "Q://KeyAntigenicSite/ViralStruct/OtherVirus/H3N2/Data/PDB/6TZB.pdb" 
    #Kvirus_chains='ABCDEF'      # 病毒部分
    #Kreceptor_chains=''    # 受体部分
    #scanner = MultiChainImmuneScanner(pdb_file, Kvirus_chains, Kreceptor_chains,'SIA')

    #SARS
    #pdb_file = "Q://KeyAntigenicSite/ViralStruct/StructurePredictESC/SARS/6ACJ.pdb" 
    #Kvirus_chains='ABC'      # 病毒部分
    #Kreceptor_chains='D'    # 受体部分
    #scanner = MultiChainImmuneScanner(pdb_file, Kvirus_chains, Kreceptor_chains)

    #H1N1pdm
    #pdb_file = "Q://KeyAntigenicSite/ViralStruct/StructurePredictESC/H1N1/4jtv.pdb" 
    #Kvirus_chains='ABCDEFGHIJKL'      # 病毒部分
    #Kreceptor_chains=''    # 受体部分
    #scanner = MultiChainImmuneScanner(pdb_file, Kvirus_chains, Kreceptor_chains,'SIA')

    #Poliovirus1
    #pdb_file = "Q://KeyAntigenicSite/ViralStruct/StructurePredictESC/Poliovirus/3J8F.pdb" 
    #Kvirus_chains='1234'      # 病毒部分
    #Kreceptor_chains='7'    # 受体部分
    #scanner = MultiChainImmuneScanner(pdb_file, Kvirus_chains, Kreceptor_chains)

    scanner = MultiChainImmuneScanner(args.pdb, args.virus, args.receptor, args.ligand)
    scanner.run_analysis()
