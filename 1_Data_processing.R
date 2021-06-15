# R script for identification of sample-specific MS features 
# written by Thomas E. Smith
# created by TES July 2020
# last updated by TES April 2021

# setting the stage
wd <- "/Users/User/Work"
setwd(wd)
library("tidyverse")
library("stats")



### Import validated, filtered Byologic peptides list and clean up
## Import data from .csv files, downloadable from GitHub (smit4227) or MassIVE
can.red <- read.csv("canonical_reduced_PGN.csv", 
                    header = T, stringsAsFactors = F)
noncan.red <- read.csv("noncanonical_reduced_PGN.csv", 
                       header = T, stringsAsFactors = F)
all.nonred <- read.csv("all_nonreduced_PGN.csv",
                       header = T, stringsAsFactors = F)


## Combine data into a single data frame, convert to tibble
all.dat <- rbind(can.red, noncan.red, all.nonred) %>%
  as_tibble()


## Remove columns that don't provide any useful information
all.dat <- all.dat %>%
  select( -c(Validate, Comment, Obs..m.z, Calc..m.z, Glycans,
             Fragment.type.s., Score, PID, Labels, MS2.Search.Alias.name, 
             Protein.isoelectric.point, Sample.type) )


## Filter out unnecessary rows
# Keepers include rows with column "Row." values containing only one decimal
all.dat <- all.dat %>%
  filter( str_count(Row., pattern = "\\.") == 1 )  %>%
  # Remove all rows containing Angiotensin II data, an internal standard that
  # was incorrectly added to certain samples
  filter( Sequence != "-.DRVYIHPF.-" )


## Define sample type
# Add a column that assigns sample type to each replicate based on alias
sample.names <- c("1A", "2A", "3A",
                  "7B", "8B", "9B",
                  "16B", "17B", "18B",
                  "blank1", "blank2", "blank3", "blank4", "blank5",
                  "sonic3", "lys3")
sample.types <- c("mut.wall", "mut.wall", "mut.wall",
                  "EcLdcA.muro", "EcLdcA.muro", "EcLdcA.muro",
                  "ApLdcA.muro", "ApLdcA.muro", "ApLdcA.muro",
                  "blank", "blank", "blank", "blank", "blank",
                  "sonication.Ap", "lysozyme.Ap")
all.dat <- all.dat %>%
  mutate( Sample.type = sample.types[ match(MS.Alias.name, sample.names) ] )


## Add a column identifying the multimeric state of the compound
# Monomer, dimer, or trimer
all.dat <- all.dat %>%
  mutate( Multimer = str_count(Sequence, pattern = "J") + str_count(Sequence, pattern = "j") )


## Add a compound description/label for each unique compound
# Create a concise description of each compound from Protein.name, Sequence,
# Glycans, and Mod...Summary columns

all.dat <- all.dat %>%
  # First, obtain stem peptide sequence from Protein.name column, listing donor
  # before acceptor (I listed acceptor before donor in Byonic/Byologic)
  mutate( Peptide = case_when( 
    # Monomers with "mDap" in place of "m"
    Multimer == 1 & str_count(Protein.name, "_") == 0 ~ 
      gsub("-mDap", "m", Protein.name),
    # Monomers with two possible sequences
    Multimer == 1 & str_count(Protein.name, "_") == 1 ~
      paste0("AEm(", 
             gsub("AEm.*_AEm", "", Protein.name), 
             "/", 
             gsub("AEm", "", gsub("_.*", "", Protein.name)),
             ")"),
    # Tetra-tetra and tri-penta dimers were not resolved
    Multimer == 2 & str_count(Protein.name, "_") == 1 &
      Protein.name == "AEmA-AEmA_AEmAA-AEm" ~ paste0("AEm(A/)-AEm(A/AA)"),
    # Tetra-tri and tri-tetra were not resolved (4-3 vs. 3-3 crosslinkage)
    Multimer == 2 & str_count(Protein.name, "_") == 1 &
      Protein.name == "AEmA-AEm_AEm-AEmA" ~ paste0("AEm(A/)-AEm(/A)"),
    # AEmA-AEmG and AEm-AEm(Q/AG) can't be resolved but one sequence was
    # used to identify some nonreduced masses while the other sequence was
    # used to identify their corresponding reduced masses. Make the 
    # nonreduced/reduced mass pair have the same peptide sequence
    Multimer == 2 & Protein.name == "AEmG-AEmA" | 
      Protein.name == "AEmAG-AEm_AEmQ-AEm" ~ paste0("AEm(A//)-AEm(G/Q/AG)"),
    # Remaining dimers listed as acceptor-donor
    Multimer == 2 & str_count(Protein.name, "_") == 0 &
      Protein.name != "AEmG-AEmA" ~ paste0("AEm", 
                                           gsub("AEm.*-AEm", "", Protein.name), 
                                           "-AEm", 
                                           gsub("AEm", "", gsub("-.*", "", Protein.name))
      ),
    # Remaining acceptor-donor dimers with two possible sequences
    Multimer == 2 & str_count(Protein.name, "_") == 1 &
      Protein.name != "AEmA-AEmA_AEmAA-AEm" &
      Protein.name != "AEmA-AEm_AEm-AEmA" &
      Protein.name != "AEmAG-AEm_AEmQ-AEm" ~
      paste0("AEm",
             gsub("AEm.*-AEm", "", Protein.name),
             "-AEm(",
             gsub("-AEm.*", "", gsub("AEm.*_AEm", "", Protein.name)),
             "/",
             gsub("AEm", "", gsub("-.*", "", Protein.name)),
             ")"),
    # Trimers listed as acceptor-middle-donor
    Multimer == 3 & str_count(Protein.name, "_") == 0 ~
      paste0("AEm",
             gsub("AEm.*-AEm", "", Protein.name),
             "-AEm",
             gsub("AEm", "", gsub("^.*-(.*)-.*", "\\1", Protein.name)),
             "-AEm",
             gsub("AEm", "", gsub("-.*", "", Protein.name))
      ),
    # Acceptor-middle-donor trimers with two possible sequences
    Multimer == 3 & str_count(Protein.name, "_") == 1 ~
      paste0("AEm(",
             gsub("AEm.*-AEm", "", gsub("_.*", "", Protein.name)),
             "/",
             gsub("AEm.*-AEm", "", gsub(".*_", "", Protein.name)),
             ")-",
             gsub("^.*-(.*)-.*", "\\1", Protein.name),
             "-AEm(",
             gsub("AEm", "", gsub("-.*", "", gsub("_.*", "", Protein.name))),
             "/",
             gsub("AEm", "", gsub("-.*", "", gsub(".*_", "", Protein.name))),
             ")") ) ) %>%
  # Globally replace all "AEm" with "AEJ", a better way to represent mDap #####
  mutate( Peptide = gsub("AEm", "AEJ", Peptide) ) %>% ######
  # Append glycans to the peptide sequence
  mutate( Compound = case_when( 
    # No glycans at all
    str_count(Mod..Summary, "480\\.1955") == 0 &
      str_count(Mod..Summary, "460\\.1693") == 0 &
      str_count(Mod..Summary, "478\\.1799") == 0 ~ paste0(Peptide),
    # One of the three possible glycans
    str_count(Mod..Summary, "480\\.1955") == 1 &
      str_count(Mod..Summary, "460\\.1693") == 0 &
      str_count(Mod..Summary, "478\\.1799") == 0 ~ paste0("GM-", Peptide),
    str_count(Mod..Summary, "480\\.1955") == 0 &
      str_count(Mod..Summary, "460\\.1693") == 0 &
      str_count(Mod..Summary, "478\\.1799") == 1 ~ paste0("GM-", Peptide),
    str_count(Mod..Summary, "460\\.1693") == 1 &
      str_count(Mod..Summary, "480\\.1955") == 0 &
      str_count(Mod..Summary, "478\\.1799") == 0 ~ paste0("GaM-", Peptide),
    # Two different glycans 
    str_count(Mod..Summary, "480\\.1955") == 1 &
      str_count(Mod..Summary, "460\\.1693") == 1 ~ paste0("GM-", Peptide, "-GaM"),
    str_count(Mod..Summary, "478\\.1799") == 1 &
      str_count(Mod..Summary, "460\\.1693") == 1 ~ paste0("GM-", Peptide, "-GaM"),
    # Two of the same glycan
    str_count(Mod..Summary, "480\\.1955") == 2 ~ paste0("GM-", Peptide, "-GM"),
    str_count(Mod..Summary, "460\\.1693") == 2 ~ paste0("GaM-", Peptide, "-GaM"),
    str_count(Mod..Summary, "478\\.1799") == 2 ~ paste0("GM-", Peptide, "-GM") )) %>%
  # Remove columns that are no longer necessary
  select( -c(Row., Protein.name, Mod..Summary) )


## Some samples were reduced and others were not. 
# Add a column identifying compounds as reduced or nonreduced by Sample.type
# and by the glycans in Compound; only compounds for nonereduced samples 
# and containing GM will be labeled as nonreduced.
nonreduced.samples <- c("sonication.Ap", "lysozyme.Ap")
all.dat <- all.dat %>%
  group_by( Calc.M ) %>%
  mutate( Compound.state = ifelse(any(Sample.type %in% nonreduced.samples) &
                                    grepl("GM", Compound),
                                  "nonreduced",
                                  "reduced") ) %>%
  ungroup()



### Adjustment for dilution
## Some samples were diluted prior to running on LC-MS. Multiply by dilution
## factor to account for this.
diluted.samples <- c("mut.wall", "EcLdcA.muro", "ApLdcA.muro")
all.dat$XIC.AUC <- as.numeric(all.dat$XIC.AUC)
all.dat <- all.dat %>%
  mutate( XIC.AUC = ifelse( Sample.type %in% diluted.samples, 
                            XIC.AUC*5, XIC.AUC) )



### Baseline subtraction
## For each compound, subtract the average AUC of the blank samples
adj.dat <- all.dat %>%
  # Calculate the average AUC of the blank samples for each compound, 
  # accounting for the different blanks for reduced and nonreduced samples 
  group_by( Calc.M, Compound.state ) %>%
  mutate( Avg.blanks = mean( XIC.AUC[grepl("blank", MS.Alias.name)] ) ) %>%
  # Subtract these values from each sample for each compound 
  mutate( Adj.AUC = XIC.AUC - Avg.blanks ) %>%
  # Change any resulting negative values to 0
  mutate( Adj.AUC = ifelse(Adj.AUC < 0, 0, Adj.AUC) ) %>%
  ungroup() %>%
  # Remove blank sample data
  filter( Sample.type != "blank" )



### Normalization
## Values should be normalized by total PGN abundance per sample to
## eliminate variation between replicates

# Calculate the sum of AUCs for each compound per sample, then divide each AUC
# value by the sample-specific sum
norm.dat <- adj.dat %>%
  group_by( MS.Alias.name ) %>%
  mutate( Total.AUC.sample = sum( Adj.AUC ) ) %>%
  mutate( Norm.AUC = Adj.AUC/Total.AUC.sample ) %>%
  ungroup()



### Separate datasets
## Separate the Ap dataset from the enzyme dataset

norm.Ap.dat <- norm.dat %>%
  filter( Sample.type == "mut.wall" | Sample.type == "sonication.Ap" | 
            Sample.type == "lysozyme.Ap" ) %>%
  mutate( Dataset = "Ap" )

norm.enzyme.dat <- norm.dat %>%
  filter( !(Sample.type == "sonication.Ap" | Sample.type == "lysozyme.Ap") ) %>%
  mutate( Dataset = "enzyme" )



### Data transformation
## We want to show the proportion of the total AUC that each individual PGN
## compound represents.
## The data are already normalized to between 0-1, so we simply remove rows
## with no data, then log transform the data

## First for the enzyme dataset
# Scale the data so that no meaningful values will be negative after 
# log transformation - these values should be greater than 1.
# To do this, first look at the distributions of the log transformed data 
# to determine an appropriate scaling factor.
hist(log10(norm.enzyme.dat$Norm.AUC), labels = T, breaks = 20, xlim = c(-9, 0))
enzyme.scaling.factor <- 1e9
# Scale, add 2, and log transform each dataset
log.enzyme.dat <- norm.enzyme.dat %>%
  mutate( Scaled.AUC = Norm.AUC * enzyme.scaling.factor ) %>%
  mutate( Scaled.AUC = Scaled.AUC +2 ) %>%
  mutate( Log.AUC = log10(Scaled.AUC) )

## Then for the Ap dataset
# Same process
hist(log10(norm.Ap.dat$Norm.AUC), labels = T, breaks = 20, xlim = c(-9, 0))
Ap.scaling.factor <- 1e9
# Scale, add 2, and log transform each dataset
log.Ap.dat <- norm.Ap.dat %>%
  mutate( Scaled.AUC = Norm.AUC * Ap.scaling.factor ) %>%
  mutate( Scaled.AUC = Scaled.AUC +2 ) %>%
  mutate( Log.AUC = log10(Scaled.AUC) )



### Save data
# export transformed datasets to .csv files
write.csv(log.enzyme.dat, 
          row.names = F,
          file = paste0("enzyme_dataset.csv"))
write.csv(log.Ap.dat, 
          row.names = F,
          file = paste0("Ap_dataset.csv"))
# save important objects as RData
save(log.enzyme.dat, log.Ap.dat,
     file = paste0("1_Data_processing.RData"))



