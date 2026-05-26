# **Structural topology of viral binding proteins dictates divergent antigenic evolutionary strategies**

## **1. Structural data collection and multidimensional scaling**

Structural data for viral binding proteins in complex with their respective receptors were retrieved from the Protein Data Bank (PDB). For the SARS-CoV-2 Spike protein, PDB IDs 7DX1, 7DX3, 7DX5, 7DX7, 7DX8, 7DX9, 7WVN, and 8HRI were utilized. For the influenza H3N2 Hemagglutinin (HA) protein, PDB IDs 2YP8, 4WEA, and 6TZB were used. Additionally, 4JTV and 6ACJ were utilized for H1N1pdm and SARS-CoV, respectively. 

Structural visualization and mapping were performed using [UCSF ChimeraX (v1.11)](https://www.cgl.ucsf.edu/chimerax/download.html). 

To represent the spatial position of each residue, the 3D coordinates of the alpha-carbon were extracted using our in-house script: __./Script/viral_structure_analyzer.py__ for each chain and then merged manually. 

For proteins comprising multiple chains or exhibiting multiple conformational states, the spatial coordinates of homologous sites were averaged across all chains and conformations prior to multidimensional scaling (MDS). The MDS spatial projection was conducted in R using the [cmdscale](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cmdscale) function.

The R script can be found in __./Script/ViralStruct.R__

## **2. Collection of functional and structural indices**

Deep mutational scanning (DMS) functional profiles were aggregated from multiple published studies (see __./Data/Supplementary Table 1.xslx__). Briefly, yeast-display measured antibody escape scores and ACE2 binding affinities were calculated from DMS data following our previously established [pipeline](https://github.com/ipplol/SARS2EVO). 

Pseudovirus-assay measured sera escape scores and Spike-mediated cell entry indices were directly obtained from datasets published by Bloom et al (The url is in the Supplementary Table 1). 

Because DMS data typically encompass multiple datasets corresponding to different viral variants, we combined them as the normalized scores of each residue:  within each dataset, we first calculated the absolute value of the mutational effect for all 20 amino acid substitutions at each site (focusing on the magnitude of the impact rather than the directional sign) and then computed the mean. Each site was subsequently normalized against the maximum value observed across the entire Spike protein. To integrate multiple variant datasets, we assigned the maximum recorded normalized score to each residue across all variants. This approach minimizes biases introduced by uneven antibody distributions across different variant exposures, such as immunological imprinting, reflecting the potential effect of each site.

Structural indices were also calculated using the in-house Python script __./Script/viral_structure_analyzer.py__ utilizing the Biopython package based on the PDB structures mentioned above. During the calculation of intrinsic structural properties, such as Weighted Contact Number (WCN) and Relative Solvent Accessibility (RSA), receptor chains were deliberately removed from the complexes. The 'distance to receptor' was defined as the minimum Euclidean distance from a given viral residue's alpha-carbon to any receptor alpha-carbon atom. Thus, users need to specify the virus chains and receptor chains.

**The final summary file including all calculated indices can be find in**
  __./Result/__ 

## **Calculation of evolutionary indices**

The mutation incidence, representing the number of independent mutational events (homoplasy), was calculated following our previously described [methodology](https://github.com/ipplol/SARS2EVO). For SARS-CoV-2, the calculation utilizes the global phylogenetic tree from the [UShER database](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/). For influenza H3N2, the calculation utilizes  a phylogenetic tree we build with [IQ-TREE](https://iqtree.github.io/) which including 5,697 complete sequences (__./Data/Supplementary Table 2.xslx__) obtained from [GISAID](https://gisaid.org/).

## **Structure-based prediction of site competitiveness**

To predict antigenic escape sites based on spatial steric hindrance, we developed a custom Python script __./Script/immune_scan.py__. The algorithm quantifies the structural "competitiveness" of each viral residue, identifying surface sites where antibody binding would physically clash with the host receptor. For each viral residue, the algorithm centers a virtual sphere with a default radius of 45Å on its alpha-carbon (which can be adjusted by user). 

In the case of receptor binding, the virtual sphere is located on the outer surface of the residue, perpendicular to the normal. While in the case of ligand binding, the virtual sphere is centered on the alpha-carbon of each residue.  It then calculates a "clash number" by counting the total number of receptor atoms located within this sphere, effectively simulating the steric bulk an antibody would encounter. To account for the baseline physical availability of the residue to antibodies, this clash number is multiplied by the residue's solvent-accessibility surface area (SASA).

**Usage:**

**immune_scan.py [-h] --pdb PDB --virus VIRUS --receptor RECEPTOR [--ligand LIGAND]**

Multi-Chain Viral Immune Escape Scanner

optional arguments:
 -h, --help        show this help message and exit
 --pdb PDB          Path to PDB file
 --virus VIRUS        Chain IDs of Virus (e.g., 'ABC' for trimer, 'A' for monomer)
 --receptor RECEPTOR  Chain IDs of Receptor (e.g., 'D', or 'HL' for antibody)
 --ligand LIGAND      [Optional] Ligand Name (e.g., SIA). If ignored, assumes Protein Receptor.

**H3N2 ligand model**
python immune_scan.py --pdb /mnt/q/KeyAntigenicSite/ViralStruct/OtherVirus/H3N2/Data/PDB/6TZB.pdb --virus ABCDEF --receptor **''** --ligand SIA

**SARS-CoV-2 protein model**
python immune_scan.py --pdb /mnt/q/KeyAntigenicSite/ViralStruct/OtherVirus/H3N2/Data/PDB/6TZB.pdb --virus ABC --receptor D


## **Machine learning prediction using EVEscape**

To benchmark our structure-based competitiveness metric against generalized sequence-based machine learning approaches, we utilized the [EVEscape model](https://evescape.org/). For the influenza H3N2 HA protein, predictions were generated locally by running the official EVEscape pipeline with default parameters, using A/Massachusetts/18/2022 as the reference sequence without any supplementary modifications or fine-tuning. For SARS-CoV-2, the EVEscape prediction scores for the Spike protein were directly retrieved from the pre-computed datasets published in the [original EVEscape study](https://github.com/OATML-Markslab/EVEscape/tree/main/results/summaries).



