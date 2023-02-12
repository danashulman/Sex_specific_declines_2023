# Sex_specific_declines_2023 
This project supports the work in 'Sex-specific declines in cholinergic-targeting tRNA fragments in the nucleus accumbens in Alzheimer’s disease' by Shulman et al. 

# Summary:
Females with Alzheimer’s disease (AD) present accelerated damage to cholinergic neurons and undergo accelerated cognitive decline compared to males, but the underlying mechanisms are incompletely understood. Here, we report that tRNA fragments (tRFs) carrying complementary sequences to cholinergic transcripts (CholinotRFs) may causally contribute to both these phenomena. We found that the nucleus accumbens of female brains, a brain region enriched in cholinergic neurons, exhibits larger CholinotRF decreases compared to hypothalamic or cortical tissues. These CholinotRFs are largely of mitochondrial genome origin, and their reduced levels are correlated with elevated levels of their predicted cholinergic associated mRNA targets. Additionally, AD temporal cortices showed sex-specific changes in distinct cell populations’ levels of cholinergic transcripts. Moreover, cholinergic differentiation of human-originated neuronal cell lines was accompanied by sex-specific elevation of CholinotRFs, supporting the role of CholinotRFs in cholinergic regulation. Our findings highlight possible CholinotRF involvement in the AD sex-specific cholinergic loss and cognitive deterioration. 

# Files:
README.md, Single cell.py, Statistical Analysis.py, Brain_Data_README.md, Cell_Lines_README.md, Single_Cell_README.md, sncRNA_Targets_README.md.

# Directories:
1. Brain/Data: Brain_Data_README.md,  normalized tRF count tables,  normalized miRNA count tables,  normalized mRNA count tables, raw tRF count tables,  raw miRNA count tables,  raw mRNA count tables,  tRF metadata from alignement to MintMap.
2. Brain/Patients_Info: Contains the metadata of the donors from ROS-MAP poject.
2. Brain/Results: FDR corrected Kolmogorov-Smirnov test results, separated by group (all/females/males), tissue (nuc/hyp/STG), type (tRF/miRNA) and trait (braak/cerad/cogdx).
4. Cell_Lines: Contains the data and results of the analysis of the human-originated cell lines. Includes the files: Cell_Lines_README.md, tRF_raw_counts.csv, tRF_meta.csv, tRF_norm_counts_2.csv, tRF_norm_counts_4.csv, info.csv, Female_time2_ttest.csv, Male_time2_ttest.csv, Female_time4_ttest.csv, Male_time4_ttest.csv.
5. Single_Cell: Single cell.py, Single_Cell_README.md, metadata.csv, n_samples.csv, percentages.csv.
6. Single_Cell/Utest: The results of the Mann-Whitney U-test on cholinergic mRNAs in each cell population.
7. Single_Cell/LFC: The LFC of each cholinergic gene as a function of linear combinations representing the MTG scRNA data.
8. tRF_Targets: Contains the targets of all the tRFs analysed in the paper.
9. miRNA_Targets: Contains the targets of all the miRs analysed in the paper.
