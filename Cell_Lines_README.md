# Sex_specific_declines_2023 / Cell-Lines
This project supports the work in 'Sex-specific declines in cholinergic-targeting tRNA fragments in the nucleus accumbens in Alzheimer’s disease' by Shulman et al.

Here, the data was used to challenge role of CholinotRFs in cholinergic regulation. we studied changes in small RNA-seq transcripts in human neuroblastoma cell lines  from  male and female origins (LA-N-2, LA-N-5) 2 and 4 days after exposure to neurokines that induce cholinergic differentiation.

# Data discription:
The Cell Lines LA-N-2 (female) (DSMZ Cat# ACC-671, RRID:CVCL_1829) and LA-N-5 (male) (DSMZ Cat# ACC-673, RRID:CVCL_0389) were purchased at DSMZ (Braunschweig, Germany). These cells respond to differentiation by several neurokines (CNTF, LIF, IL-6) by cholinergic differentiation according to a known protocol97, corresponding with elevation of choline acetyltransferase (ChAT), the central cholinergic marker (mRNA, protein, and activity) as well as it’s intronic vesicular ACh-transporter gene, vAChT (aka SLC18A3). The small-RNA sequencing files were produced in-house and are available to all interested researchers from the GEO data portal (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132951, GEO accession: GSE13295). Each time point (2 days vs 4 days) contained 4 replicates from each state (Control vs Differentiated) for each sex (Female vs Male). 

# Files:
1. tRF_raw_counts.csv: The raw counts of the tRFs aligned by MintMap, downloadd from NIH.
2. tRF_meta.csv: The metadata of the tRFs after the alignment.
3. tRF_norm_counts_2.csv: tRF counts 2 days after exposure, normalized using Deseq2 'median of ratios' method
4. tRF_norm_counts_4.csv: tRF counts 4 days after exposure, normalized using Deseq2 'median of ratios' method.
5. info.csv: The metadata of the cell lines.
6. Female_time2_ttest.csv: t-test results in female-originated cell lines 2 days after exposure. 
7. Male_time2_ttest.csv: t-test results in male-originated cell lines 2 days after exposure. 
8. Female_time4_ttest.csv: t-test results in female-originated cell lines 4 days after exposure. 
9. Male_time4_ttest.csv: t-test results in male-originated cell lines 4 days after exposure. 

The analyses were prformed using the code in 'Statistical Analysis.py'.
