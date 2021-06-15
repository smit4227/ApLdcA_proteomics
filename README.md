# ApLdcA_proteomics

Scripts for processing and analysis of 1) Buchnera muropeptides isolated from aphids, and 2) E. coli muropeptides treated with the A. pisum LdcA enzyme, ApLdcA

Thomas E. Smith, University of Texas at Austin, Department of Intregrative Biology, Nancy A. Moran lab

These scripts are associated with the publication, "...", doi:... . Raw LC-MS/MS data is available from the MassIVE database (...). In this study, we characterized composition of peptidoglycan muropeptides isolated from aphids and derived from the aphid symbiont Buchnera aphidicola. Additionally, we produced and characterized the aphid enzyme ApLdcA by treating muropeptides derived from E. coli cell walls with ApLdcA or the E. coli homolog, EcLdcA. 

What you need to repeat our analysis:

The "Reference Files" folder contains the output files from Byologic, the proteomics software used to identify and quantify muropeptides. The three files are named based on whether they contain data for muropeptides with canonical (e.g. AE-mDap-AA) or noncanonical (e.g. AE-mDap-G) stem peptides, and whether or not the anomeric centers of the muropeptide carbohydrates were reduced with sodium borohydride or not (aphid muropeptides were not reduced, while all E. coli muropeptides were reduced). We used the following programs in our analyses:
R programs: R v3.6.1 RStudio v1.1.463 tidyverse v1.3.0 stats v3.6.1 broom v0.7.2 colorspace v1.4-1 ggplot2 v3.2.1 ggpubr v0.4.0

What to do:

Download the reference files and "1_Data_processing.R" into the same directory and run the script to generate a curated spreadsheet of muropeptides and their relative abundance per sample. These output files are used for the input of "2_Stats_and_Plots.r", the output of which are the plot files used to construct Figures 2-3, S1, S3 in the paper.

Updated June 2021
