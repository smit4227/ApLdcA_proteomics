# ApLdcA_proteomics

Scripts for processing and analysis of 1) Buchnera muropeptides isolated from aphids, and 2) E. coli muropeptides treated with the A. pisum LdcA enzyme, ApLdcA

Thomas E. Smith, University of Texas at Austin, Department of Intregrative Biology, Nancy A. Moran lab

These scripts are associated with the publication, "...", doi:... . Raw LC-MS/MS data is available from the MassIVE database (...). In this study, we characterized composition of peptidoglycan muropeptides isolated from aphids and derived from the aphid symbiont Buchnera aphidicola. Additionally, we produced and characterized the aphid enzyme ApLdcA by treating muropeptides derived from E. coli cell walls with ApLdcA or the E. coli homolog, EcLdcA. 

What you need to repeat our analysis:

The "Reference Files" folder contains the output files from Byologic, the proteomics software used to identify and quantify muropeptides. The three files are named based on whether they contain data for muropeptides with canonical (e.g. AE-mDap-AA) or noncanonical (e.g. AE-mDap-G) stem peptides, and whether or not the anomeric centers of the muropeptide carbohydrates were reduced with sodium borohydride or not (aphid muropeptides were not reduced, while all E. coli muropeptides were reduced). We used the following programs in our analyses:
R programs: R v3.6.1 RStudio v1.1.463 tidyverse v1.3.0 stats v3.6.1 broom v0.7.2 colorspace v1.4-1 ggplot2 v3.2.1 ggpubr v0.4.0

What to do:

Download the reference files and "1_Data_processing.R" into the same directory and run the script to generate a curated spreadsheet of muropeptides and their relative abundance per sample. These output files are used for the input of "2_Stats_and_Plots.r", the output of which are the plot files used to construct Figures 2-3, S1A-B,E-G in the paper, for which the data are contained in the .xslx file described below.

Legends for tabs A-C of the file "Muropeptide_raw_transformed_stat_data.xlsx":

Quantification data (A and B) and statistical analyses (C) from the proteomic-analysis of A) muropeptides isolated from A. pisum and E. coli (Figures 2, S1, Table S1) and B-C) LdcA-treated E. coli muropeptides (Figures 3, S3, Table S1). A-B) Samples are indicated by the “MS.Alias.name” column (sonic3 and lys3 for A. pisum; 1A, 2A, and 3A for untreated E. coli, 7B, 8B, and 9B for EcLdcA-treated E. coli, and 16B, 17B, 18B for ApLdcA-treated E. coli samples). The “Sample.type” column describes the origin of the sample in the format of enzyme/treatment (sonication or lysozyme treatment for A. pisum; mut=mutanolysin for E. coli) and substrate (Ap for A. pisum muropeptides; wall for E. coli sacculus). Stem peptide sequences are indicated by the “Sequence” and “Peptide” columns, where “J” represents mDap. The “Compound” adds glycan substituents to the stem peptide sequence. The “Multimer” column indicates the number of stem peptides per compound. “Compound.state” refers to whether or not the sample was reduced with sodium borohydride prior to HPLC purification. “Adj.AUC” is equivalent to “XIC.AUC” minus “Avg.blanks” within each row. For each row, “Norm.AUC” is the proportion of total PGN that the corresponding compound represents, and is equivalent to “Adj.AUC” divided by “Total.AUC.sample” (averaged values per sample type shown in Table S1). Data were then scaled by adding 1e9 to each “Norm.AUC” value (“Scaled.AUC”) to avoid negative numbers during log transformation, and 2 was added to each value to avoid non-real numbers during log transformation. The final “Log.AUC” data are shown in Figures 2-3, S1, S3. C) Tukey’s Honest Significant Difference (HSD) test results of the LdcA-treated E. coli muropeptide dataset (Figure 3). For each of 90 PGN compounds, the HSD test was applied for each sample type pair (three pair-wise comparisons per compound). The resulting p-values were adjusted for false-discovery rate via multiple test correction, and statistical significance attributed to comparisons with “fdr.p.value” < 0.05.

Updated July 2021
