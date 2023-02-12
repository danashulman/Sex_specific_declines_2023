# Sex_specific_declines_2023 / Single_Cell

This project supports the work in 'Sex-specific declines in cholinergic-targeting tRNA fragments in the nucleus accumbens in Alzheimer’s disease' by Shulman et al.

Here, the data was used to identify the cell types that contributed to the altered levels of cholinergic mRNAs seen in these extra NAc regions. Specifically, we used single-cell RNA-Sequencing (scRNA-Seq) data from the Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) study, to determine cell type-specific changes in cortical neurons. 

# Data discription:
The scRNA analysis was performed on cellular level transcriptomic data from the middle temporal gyrus (MTG) of female and male aged volunteers on the AD spectrum. The data were downloaded from the Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD), and can be accessed from https://registry.opendata.aws/allen-sea-ad-atlas. Study data were generated from postmortem brain tissue obtained from the MTG of 84 aged individuals spanning the full spectrum of AD severity (preMCI, MCI, mild-moderate severe AD) and 5 neurotypical aged adult cognitively intact individuals. The data included classification to 24 distinct cell populations, of which 18 were of neuronal origin and 6 of non-neuronal origin (astrocytes, oligodendrocytes, OPCs, microglia-PVM, endothelial and VLMC). Notably, the non-neuronal cells constituted only ~20% (20.7% in females, 19.3% in males) of all the cells in the analysis.  

# Files:
1. Single cell.py: The code for the analyses performed on the MTG scRNA data.
2. metadata.csv: The metadata of the donors.
3. n_samples.csv: The number of cells in each analysis for each specific cell population.
4. percentages.csv: The relative proportion of each cell population from all the cells.
5. Single_Cell_Utest directory: The results of the Mann-Whitney U-test on cholinergic mRNAs in each cell population.
6. Single_Cell_LFC directory: The LFC of each cholinergic gene as a function of linear combinations representing the MTG scRNA data.
