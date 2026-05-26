import numpy as np
import pandas as pd
import re
import warnings
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBParser, Select, PDBIO
from Bio.PDB.Polypeptide import is_aa, three_to_one
from scipy.spatial import distance_matrix

# 忽略 Biopython 中不影响计算的警告
warnings.filterwarnings('ignore')

class ViralMultimerAnalyzer:
    def __init__(self, pdb_file, virus_chains, receptor_chains=None):
        #"""
        #初始化分析器，明确区分病毒本体和受体
        #
        #:param pdb_file: PDB 文件路径
        #:param virus_chains: list, 病毒的多聚体链 ID 列表 (e.g., ['A', 'B', 'C'])
        #:param receptor_chains: list, 受体的链 ID 列表 (e.g., ['D', 'E']), 如果没有则为 None
        #"""
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("complex", pdb_file)
        self.model = self.structure[0]
        self.pdb_file = pdb_file
        
        self.v_chains = virus_chains
        self.r_chains = receptor_chains if receptor_chains else []
        
        print(f"Loaded {pdb_file}")
        print(f"  - Virus Chains (Trimer/Oligomer): {self.v_chains}")
        print(f"  - Receptor Chains (Excluded from WCN/RSA): {self.r_chains}")

    def _get_coords(self, chain_ids, atom_type='CA'):
        #"""辅助函数：获取指定链列表中所有指定原子的坐标"""
        coords = []
        for cid in chain_ids:
            if cid in self.model:
                chain = self.model[cid]
                for res in chain:
                    if is_aa(res) and atom_type in res:
                        coords.append(res[atom_type].get_coord())
        return np.array(coords)

    def calculate_rsa(self):
        #"""
        #计算 RSA (相对溶剂可及性)。
        #【关键逻辑】：为了获得病毒表面真实的抗体可及性，我们需要在一个
        #“移除受体”的临时模型上运行 DSSP。否则受体覆盖的区域会被误判为掩埋。
        #"""
        # 1. 创建一个临时的 PDB 文件，只包含病毒链
        class VirusSelect(Select):
            def __init__(self, v_chains): self.v_chains = v_chains
            def accept_chain(self, chain): return chain.get_id() in self.v_chains
            
        io = PDBIO()
        io.set_structure(self.model)
        temp_file = "temp_virus_only.pdb"
        io.save(temp_file, VirusSelect(self.v_chains))
        
        # 2. 在只有病毒的结构上运行 DSSP
        # 注意：这里 DSSP 会把 A+B+C 作为一个整体，正确处理三聚体界面
        try:
            p = PDBParser(QUIET=True)
            temp_struct = p.get_structure("temp", temp_file)
            dssp = DSSP(temp_struct[0], temp_file, dssp='mkdssp')
            
            rsa_dict = {}
            for key in dssp.keys():
                chain_id = key[0]
                res_id = key[1][1]
                # key[3] 是 RSA
                rsa_dict[(chain_id, res_id)] = dssp[key][3]
            return rsa_dict
        except Exception as e:
            print(f"Error calculating RSA (check DSSP installation): {e}")
            return {}

    def calculate_wcn(self, target_chain_id):
       # """
       # 计算 WCN (加权接触数)。
       # 【关键逻辑】：
       # - Target: 目标链的残基 (e.g., Chain A)
       # - Environment: 整个病毒多聚体 (Chain A + B + C)
       # - Exclude: 受体链 (Receptor) 不参与计算，因为我们要看的是病毒自身的结构约束。
       # """
        # 1. 获取环境坐标 (所有病毒链的 CA 原子)
        env_coords = self._get_coords(self.v_chains, 'CA')
        
        # 2. 获取目标链残基
        target_chain = self.model[target_chain_id]
        target_residues = [r for r in target_chain if is_aa(r) and 'CA' in r]
        target_coords = np.array([r['CA'].get_coord() for r in target_residues])
        target_ids = [r.get_id()[1] for r in target_residues]
        
        if len(env_coords) == 0 or len(target_coords) == 0:
            return {}

        # 3. 计算距离矩阵 (Target vs All Virus)
        d_mat = distance_matrix(target_coords, env_coords)
        
        # 4. 自身距离设为无穷大
        d_mat[d_mat == 0] = np.inf
        
        # 5. 计算 WCN = sum(1/d^2)
        wcn_values = np.sum(1 / (d_mat ** 2), axis=1)
        
        return dict(zip(target_ids, wcn_values))

    def calculate_dist_to_receptor(self, target_chain_id):
        #"""
        #计算到受体的最小距离。
        #只有在这里，受体链才会被用到。
        #"""
        if not self.r_chains:
            return {}
            
        # 1. 获取受体所有原子坐标 (作为“墙”)
        # 这里用所有原子而不仅是CA，为了更精确的物理接触判断
        receptor_atoms = []
        for cid in self.r_chains:
            if cid in self.model:
                for atom in self.model[cid].get_atoms():
                    receptor_atoms.append(atom.get_coord())
        
        if not receptor_atoms: return {}
        receptor_coords = np.array(receptor_atoms)
        
        # 2. 遍历目标链残基
        dist_dict = {}
        target_chain = self.model[target_chain_id]
        
        for res in target_chain:
            if not is_aa(res): continue
            
            # 提取该残基所有原子
            res_atoms = np.array([a.get_coord() for a in res.get_atoms()])
            
            # 计算该残基任意原子到受体任意原子的最小距离
            dists = distance_matrix(res_atoms, receptor_coords)
            min_dist = np.min(dists)
            dist_dict[res.get_id()[1]] = min_dist
            
        return dist_dict

    def calculate_glyco_shield(self, target_chain_id):
        #"""
        #计算糖基化屏蔽效应。
        #【关键逻辑】：
        #糖链可能来自同一条链，也可能来自邻近的另一条病毒链 (跨链屏蔽)。
        #因此，我们要先找到所有病毒链上的糖基化位点，再计算目标链残基到它们的距离。
        #"""
        # 1. 寻找所有病毒链上的 N-linked 糖基化位点坐标
        all_glyco_coords = []
        
        for cid in self.v_chains:
            chain = self.model[cid]
            seq = ""
            res_obj_list = []
            
            # 构建序列
            for res in chain:
                if is_aa(res):
                    try:
                        seq += three_to_one(res.get_resname())
                        res_obj_list.append(res)
                    except: pass
            
            # 正则匹配 N-X-S/T
            # (?=(...)) 用于捕获重叠匹配
            for match in re.finditer(r'(?=(N[^P][ST]))', seq):
                idx = match.start()
                res = res_obj_list[idx]
                # 优先用侧链原子代表糖基挂载点
                if 'ND2' in res:
                    all_glyco_coords.append(res['ND2'].get_coord())
                elif 'CA' in res:
                    all_glyco_coords.append(res['CA'].get_coord())
                    
        if not all_glyco_coords: return {}
        all_glyco_coords = np.array(all_glyco_coords)
        
        # 2. 计算目标链每个残基到最近糖基化点的距离
        dist_dict = {}
        target_chain = self.model[target_chain_id]
        
        for res in target_chain:
            if is_aa(res) and 'CA' in res:
                res_coord = res['CA'].get_coord()
                # 欧氏距离
                dists = np.linalg.norm(all_glyco_coords - res_coord, axis=1)
                dist_dict[res.get_id()[1]] = np.min(dists)
                
        return dist_dict

    def calculate_bfactor_z(self, target_chain_id):
        #"""
        #计算归一化的 B-factor (Z-score)。
        #仅基于目标链自身的统计分布。
        #"""
        chain = self.model[target_chain_id]
        bfactors = []
        ids = []
        
        for res in chain:
            if is_aa(res) and 'CA' in res:
                bfactors.append(res['CA'].get_bfactor())
                ids.append(res.get_id()[1])
        
        if not bfactors: return {}
        
        vals = np.array(bfactors)
        mean = np.mean(vals)
        std = np.std(vals)
        
        if std == 0: return dict(zip(ids, np.zeros_like(vals)))
        
        z_scores = (vals - mean) / std
        return dict(zip(ids, z_scores))

    def run_full_analysis(self, target_chain_id):
        #"""
        #执行所有分析并返回 DataFrame
        #"""
        print(f"Analyzing Target Chain: {target_chain_id}...")
        
        # 1. 运行各项计算
        rsa_data = self.calculate_rsa() # 这是一个(chain, id)的字典
        wcn_data = self.calculate_wcn(target_chain_id)
        recep_dist_data = self.calculate_dist_to_receptor(target_chain_id)
        glyco_dist_data = self.calculate_glyco_shield(target_chain_id)
        bfactor_data = self.calculate_bfactor_z(target_chain_id)

        # 2. 组装结果
        results = []
        chain = self.model[target_chain_id]
        
        for res in chain:
            if not is_aa(res): continue
            rid = res.get_id()[1]
            rname = res.get_resname()
            
            row = {
                "Chain": target_chain_id,
                "Residue_ID": rid,
                "AA": rname,
                #"RSA_Unbound": rsa_data.get((target_chain_id, rid), np.nan),
                "WCN_VirusOnly": wcn_data.get(rid, np.nan),
                "Dist_to_Receptor": recep_dist_data.get(rid, np.nan),
                "Dist_to_Glycan": glyco_dist_data.get(rid, np.nan),
                "B_factor_Z": bfactor_data.get(rid, np.nan)
            }
            results.append(row)
            
        return pd.DataFrame(results)

    def sigmoid(x):
        return 1 / (1 + np.exp(-x))

    def normalization(X):
        return (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))

# --- 使用示例 ---
if __name__ == "__main__":
    # 假设你有一个包含 Spike 三聚体 (A,B,C) 和 ACE2 (D) 的文件
    # PDB File needs to be downloaded first
    pdb_file = "Q://KeyAntigenicSite/ViralStruct/OtherVirus/Poliovirus1/Data/PDB/1HXS.pdb" 
    
    # 1. 实例化：明确定义谁是病毒，谁是受体
    # 注意：你需要根据你的PDB文件实际情况修改这里的链ID
    # 例如：SARS-CoV-2 Spike通常是三聚体，但有些PDB只解析了一个RBD(E链)结合ACE2(A链)
    # 假设情况 A: 完整 Spike (A,B,C) + ACE2 (D,E,F)
    # analyzer = ViralMultimerAnalyzer(pdb_file, virus_chains=['A','B','C'], receptor_chains=['D','E','F'])
    
    # 假设情况 B: 单个 RBD (E) + ACE2 (A) (这是 6M0J 的情况)
    # 虽然是单体，但逻辑依然适用
    Kvirus_chains=['1','2','3','4']      # 病毒部分
    Kreceptor_chains=['7']    # 受体部分

    try:
        analyzer = ViralMultimerAnalyzer(
            pdb_file, 
            virus_chains=Kvirus_chains,      # 病毒部分
            receptor_chains=Kreceptor_chains    # 受体部分
        )
        
        # 2. 分析每条病毒链
        for Ichain in Kvirus_chains:
            df = analyzer.run_full_analysis(target_chain_id=Ichain)
            # 3. 保存
            df.to_csv(pdb_file+"."+Ichain+".viral_features", index=False)
            print("\nPreview of the data:")
            print(df.head())
        
    except FileNotFoundError:
        print("Can find .pdb file")