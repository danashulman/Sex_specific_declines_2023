# Sex_specific_declines_2023 / miRNA_Targets
This project supports the work in 'Sex-specific declines in cholinergic-targeting tRNA fragments in the nucleus accumbens in Alzheimer’s disease' by Shulman et al.
Predicted targets of the miRNAs by the DIANA microT prediction tool.

# Data discription:
The predicted targets of the miRs were found using the DIANA-microT-CDS prediction algorithm, to enable the best comparison between tRF and miR targets. Predicted targets with a prediction score < 0.8 were extracted. The python libraries selenium and multiprocessing were used for web scraping in order to automate the use of the DIANA prediction tools.

# Directory:
miRNA_Targets: Directory that contains the targets of all the miRs analysed in the paper.
