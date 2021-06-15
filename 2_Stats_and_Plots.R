# R script for identification of sample-specific MS features 
# written by Thomas E. Smith
# created by TES July 2020
# last updated by TES April 2021

# setting the stage
wd <- "/Users/User/Work"
setwd(wd)
library("tidyverse")
library("broom")
library("stats")
library("colorspace")
library("ggplot2")
library("ggpubr")



### Import processed LC-MS/MS data...
## ...generated from running "1_Data_processing.R" script
lnames <- load("1_Data_processing.RData")
lnames
## ... or by loading from Datasets 1 and 2 (csv files)
log.enzyme.dat <- read.csv("enzyme_dataset.csv", 
                           header = T, stringsAsFactors = F)
log.Ap.dat <- read.csv("Ap_dataset.csv", 
                       header = T, stringsAsFactors = F)



### Setup some generic plotting parameters

# Define categories and plotting order of PGN compounds, with compounds 
# containing GM coming before those with GaM, and GaM before those with none,
# and compounds containing two glycans always before those with one, and one 
# before none.
can.mono <- c("GM-AEJ", "GaM-AEJ", 
              "GM-AEJA", "GaM-AEJA", 
              "GM-AEJAA", "GaM-AEJAA")
ncan.mono <- c("GM-AEJG", "GaM-AEJG", 
               "GM-AEJS", 
               "GM-AEJV",  
               "GM-AEJC", 
               "GM-AEJ(I/L)",
               "GM-AEJ(N/GG)", "GaM-AEJ(N/GG)", 
               "GM-AEJ(Q/AG)", 
               "GM-AEJK", "GaM-AEJK", 
               "GM-AEJE", 
               "GM-AEJM", 
               "GM-AEJH", "GaM-AEJH", 
               "GM-AEJF", "GaM-AEJF", 
               "GM-AEJY", 
               "GM-AEJW", 
               "GM-AEJKR", "GaM-AEJKR")
can.di <- c("GM-AEJ-AEJ-GM", "GM-AEJ-AEJ-GaM", "GaM-AEJ-AEJ-GaM", "GM-AEJ-AEJ", 
            "GM-AEJ(A/)-AEJ(/A)-GM", "GM-AEJ(A/)-AEJ(/A)-GaM", "GaM-AEJ(A/)-AEJ(/A)-GaM", "GM-AEJ(A/)-AEJ(/A)", "GaM-AEJ(A/)-AEJ(/A)", 
            "GM-AEJ(A/)-AEJ(A/AA)-GM", "GM-AEJ(A/)-AEJ(A/AA)-GaM", "GaM-AEJ(A/)-AEJ(A/AA)-GaM", "GM-AEJ(A/)-AEJ(A/AA)", "GaM-AEJ(A/)-AEJ(A/AA)", 
            "GM-AEJA-AEJAA-GM", "GM-AEJA-AEJAA-GaM", "GM-AEJA-AEJAA", "AEJA-AEJAA") 
ncan.di <- c("GM-AEJ-AEJG-GM", "GM-AEJ-AEJG-GaM", "GM-AEJ-AEJG", 
             "GM-AEJ(A//)-AEJ(G/Q/AG)-GM", "GM-AEJ(A//)-AEJ(G/Q/AG)-GaM", "GaM-AEJ(A//)-AEJ(G/Q/AG)-GaM", "GM-AEJ(A//)-AEJ(G/Q/AG)", 
             "GM-AEJ-AEJV-GM", "GM-AEJ-AEJV-GaM", 
             "GM-AEJ-AEJ(N/GG)-GM", "GM-AEJ-AEJ(N/GG)-GaM", 
             "GM-AEJA-AEJ(N/GG)-GM", "GM-AEJA-AEJ(N/GG)-GaM", "GM-AEJA-AEJ(N/GG)", 
             "GM-AEJ-AEJD-GaM", "GaM-AEJ-AEJD-GaM", "GaM-AEJ-AEJD", 
             "GM-AEJA-AEJD-GaM", "GaM-AEJA-AEJD-GaM", "GaM-AEJA-AEJD", 
             "GM-AEJA-AEJ(Q/AG)-GM", "GM-AEJA-AEJ(Q/AG)-GaM", "GM-AEJA-AEJ(Q/AG)", 
             "GM-AEJ-AEJE",
             "GM-AEJA-AEJE", 
             "GM-AEJ-AEJF-GM", "GM-AEJ-AEJF-GaM",
             "GM-AEJA-AEJF-GM", "GM-AEJA-AEJF-GaM", "GM-AEJA-AEJF",
             "GM-AEJ-AEJKR-GM", "GM-AEJ-AEJKR-GaM", "GaM-AEJ-AEJKR-GaM", "GM-AEJ-AEJKR", "GaM-AEJ-AEJKR",
             "GM-AEJA-AEJKR-GM", "GM-AEJA-AEJKR-GaM", "GaM-AEJA-AEJKR-GaM", "GM-AEJA-AEJKR", "GaM-AEJA-AEJKR")
can.tri <- c("GM-AEJ-AEJA-AEJA-GM", 
             "GM-AEJ(A/)-AEJA-AEJ(A/AA)-GM", "GM-AEJ(A/)-AEJA-AEJ(A/AA)", "AEJ(A/)-AEJA-AEJ(A/AA)", 
             "GM-AEJA-AEJA-AEJAA", "AEJA-AEJA-AEJAA")
compound.order <- c(can.mono, ncan.mono, can.di, ncan.di, can.tri)

## Enzyme and substrate will be combined into a single parameter for plotting,
## and will be distinguished on the plot by color.
# Define an order and colors for enzyme.substrate that determines the plotting
# order of samples for each compound. Last two samples/colors are for
# background colors
sample.order <- c("mut.wall", 
                  "EcLdcA.muro", "ApLdcA.muro",
                  "A", "B")
# Colors derived from viridisLite package
sample.colors <- c("#440154FF", 
                   "#38598CFF", "#51C56AFF", 
                   NA, "gray90")
names(sample.colors) <- c(1:length(sample.order))



### Add any missing compounds to each dataset

# Add compounds missing from each dataset as a single row with 
# Sample.type = "mut.wall" and MS.Alias.name = "1A" because 
# this sample is in each dataset. Log.AUC data will be added in later.
log.enzyme.dat <- log.enzyme.dat %>%
  complete( Compound, 
            fill = list( Sample.type = "mut.wall",
                         MS.Alias.name = "1A" ) ) %>%
  # Fill in the Peptide column of the added rows with the appropriate value
  group_by( Compound ) %>%
  mutate( Peptide = first(na.omit(Peptide)) ) %>%
  ungroup()



### Add any missing sample data to the dataset  
# Each compound/sample pair should have plot points: the dataset should   
# contain X unique(compounds) * Y unique(sample.types) * N replicates rows
log.enzyme.dat <- log.enzyme.dat %>%
  complete( Compound, MS.Alias.name, 
            fill = list( Log.AUC = 0.30103 ) ) %>%
  group_by( Compound, MS.Alias.name ) %>%
  mutate( Sample.type = case_when( MS.Alias.name %in% c("1A", "2A", "3A") ~ "mut.wall",
                                   MS.Alias.name %in% c("7B", "8B", "9B") ~ "EcLdcA.muro",
                                   MS.Alias.name %in% c("16B", "17B", "18B") ~ "ApLdcA.muro"),
          Peptide = log.enzyme.dat$Peptide[match(Compound, log.enzyme.dat$Compound)] ) %>%
  ungroup()



### Convert key columns into a numeric convention

## Converting Compound, Sample.type, and Peptide columns into numbers
## will facilitate dodging at each position on the x-axis so that data
## from 3-5 samples can be plotted at each PGN compound/x-axis position. 
## Dodging resets the column order alphabetically, so converting to numbers 
## first will maintain the current order.
log.enzyme.dat <- log.enzyme.dat %>%
  mutate( Compound = as.numeric(factor(Compound, levels = compound.order)),
          Sample.type = as.numeric(factor(Sample.type, levels = sample.order)) ) %>%
  # Arrange the dataset first by compound and then by sample type
  arrange( Compound, Sample.type )
# To maintain Peptide order and facilitate determining the plot background 
# color for each peptide series, extract the Peptide order from the 
# current, sorted data frame
peptide.order <- levels(factor(log.enzyme.dat$Peptide , 
                               levels = unique(log.enzyme.dat$Peptide)))
# Convert Peptide into numbers to maintain Peptide order and facilitate
# determining the background color for each peptide series
log.enzyme.dat <- log.enzyme.dat %>%
  mutate( Peptide = as.numeric(factor(Peptide, levels = peptide.order)) )



### Perform Tukey's HSD test on each sample type pair for each compound
log.enzyme.thsd <- log.enzyme.dat %>%
  group_by( Compound ) %>% 
  do(multitst = tidy( TukeyHSD( aov(Log.AUC ~ as.factor(Sample.type),
                                    data = .) ) ) ) %>%
  unnest( multitst ) %>%
  ungroup()
# Multiple test correction by BH-FDR with only the comparisons of interest
log.enzyme.thsd <- log.enzyme.thsd %>%
  mutate( fdr.p.value = p.adjust(adj.p.value, method = "BH" ) )



### Generate plots for muropeptides, mutanolysin, and cell wall datasets
## Plot PGN compounds on x-axis, log.AUC/proportion on y-axis, and
## colors specifiy enzyme.substrate combinations.

## Split the dataset up further into typical monomers & dimers, 
# atypical monomers, atypical dimers, and trimers, to be plotted separately
df.can <- log.enzyme.dat %>%
  filter( Compound %in% match(c(can.mono, can.di), compound.order) )
df.ncan.mono <- log.enzyme.dat %>%
  filter( Compound %in% match(ncan.mono, compound.order) )
df.ncan.di <- log.enzyme.dat %>%
  filter( Compound %in% match(ncan.di, compound.order) )
df.can.tri <- log.enzyme.dat %>%
  filter( Compound %in% match(can.tri, compound.order) )

## Loop through the subsetted data frames
dfs <- list(df.can, df.ncan.mono, df.ncan.di, df.can.tri)
names(dfs) <- c("canonical_monomers_and_dimers", "noncanonical_monomers",
                "noncanonical_dimers", "canonical_trimers")

for (i in 1:length(dfs)) {
  df <- dfs[[i]]
  df.name <- names(dfs)[i]
  
  ## Establish some dataset-specific plotting parameters
  # Create custom y-axis labels for the log scale values
  y.scale <- c( floor(min(log.enzyme.dat$Log.AUC)):ceiling(max(log.enzyme.dat$Log.AUC)) )
  y.breaks <- c(0: max(y.scale))
  y.labels <- c("ND", "1E-6%", "1E-5%", "1E-4%", "0.001%", 
                "0.01%", "0.1%", "1%", "10%", "100%")
  # Define cartesian limits for y-axis, adding 2 to max y-value to make room
  # for plotting significance bars, some of which exist above the y range
  y.lims <- c(0, max(y.scale) + 2)
  # Define x-axis labels
  x.lims <- c(0, length(unique(df$Compound)))
  x.labels <- compound.order[sort(unique(df$Compound))]
  # Define colors
  color.values <- sample.colors[sort(unique(df$Sample.type))]
  color.labels <- sample.order[sort(unique(df$Sample.type))]
  
  # Establish the width of each boxplot. The actual width of each dodged 
  # boxplot will be this value divided by the number of samples plotted
  # at that point.
  box.width <- 0.75
  
  ## Differentiate between different series of compounds that all share the
  ## same stem peptide by using alternating white and gray backgrounds.
  ## This can be accomplished by plotting gray rectangles behind the dot plot
  ## at the appropriate coordinates.
  # Define the plotting order based on the Compound order of this dataset
  df <- df %>%
    mutate( plotting.order = as.numeric(factor(Compound , 
                                               levels = unique(Compound))) )
  # Determine coordinates and background colors for each peptide series
  bg.pos <- df %>%
    select( Peptide, plotting.order ) %>%
    group_by( Peptide ) %>%
    summarise( xStart = min(plotting.order), xEnd = max(plotting.order) ) %>%
    # Add a column with the peptide series number for this specific subset
    mutate( series.number = row_number() ) %>%
    # Adjust values so that borders occur between x-values
    mutate( xStart = ifelse(xStart == 1, xStart - 1, xStart - 0.5),
            xEnd = ifelse(xEnd == max(xEnd), Inf, xEnd + 0.5) ) %>%
    # Retrieve y values that have already defined
    mutate ( yStart = -Inf, yEnd = max(y.breaks) ) %>%
    # Specify whether this range is white or gray with "A and "B, respectively
    mutate( Background = ifelse( series.number %% 2 != 0, "A", "B" ) ) %>%
    ungroup()
  
  ## Generate plots
  # Calc.M on x-axis, Log.AUC on y-axis, Sample.type is color, and
  # Compound is x-axis label
  q <- ggplot() +
    # First plot the background sections as rectangles
    geom_rect(data = bg.pos, show.legend = F, 
              aes(xmin = xStart, xmax = xEnd, ymin = yStart, ymax = yEnd, 
                  fill = Background)) +
    # Then plot the data points over top of the background sections
    # Must convert Calc.M and Sample.type to factors here for some reason...
    geom_point(data = df, show.legend = F,
               aes(x = as.factor(Compound), y = Log.AUC, 
                   fill = as.factor(Sample.type)),
               size = 5, pch = 21, color = "black", 
               position = position_dodge(width = box.width)) +
    # Set the plot dimensions and point colors
    coord_cartesian(ylim = y.lims) +
    # Add colors manually, last two values are for background rectangles
    scale_fill_manual(values = c(color.values, NA, "gray90"),
                      labels = c(color.labels, NA, NA) ) +
    scale_x_discrete(labels = x.labels) +
    scale_y_continuous(breaks = y.breaks,
                       labels = y.labels) +
    # Make it pretty with some custom settings
    annotation_logticks(sides = "l", scaled = T, color = "black") +
    theme(plot.title = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 1, colour = "black"),
          axis.ticks.length = unit(7, "pt"),
          axis.line.x = element_line(color = "black", size = 1),
          axis.text = element_text(color= "black", size = 18),
          axis.text.x = element_text(size = 16, angle = 45, 
                                     hjust = 1, vjust = 1),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 0, 0, 80), "pt") )
  
  ## Add signifiance bars to the plot
  # Set some generic parameters for this section
  p.value <- 0.05
  dodge.width <- box.width/length(unique(log.enzyme.dat$Sample.type))
  comp.order <- c("1-2", "2-1", "2-3", "3-2", "1-3", "3-1")
  
  # Obtain max Log.AUC values for each Calc.M group to help determine
  # y-coordinate values later, and carry over plotting.order
  max.y <- df %>%
    group_by( Compound, plotting.order ) %>%
    summarise( Max = max(Log.AUC) ) %>%
    ungroup()
  
  # Bars must be added in layers so that they do not overlap.
  df.thsd <- log.enzyme.thsd %>%
    # Filter by the compounds contained in the current data subset
    filter( Compound %in% unique(df$Compound) ) %>%
    # Add max y values for each Compound group
    left_join( max.y, by = "Compound" ) %>%
    # Filter by significant test results, add a column for the significance
    # symbol to be plotted
    filter( fdr.p.value < p.value ) %>%
    mutate( symbol = "*" ) %>%
    # Add y-coordinates based on max Log.AUC value for each Calc.M group and
    # from the samples being compared.
    # Significance must also be layered - three possible layers
    group_by( Compound ) %>%
    mutate( y = ifelse(contrast %in% comp.order[1:4], Max + 1.3, 
                       ifelse(contrast %in% comp.order[5:6], Max + 2.0, 
                              Max + 3)) ) %>%
    # Add x-coordinates based on dataset.name, plotting.order, and contrast
    mutate( xStart = ifelse(contrast %in% comp.order[1:2], 
                            plotting.order - dodge.width,
                            ifelse(contrast %in% comp.order[3:4], 
                                   plotting.order,
                                   plotting.order - dodge.width) ),
            xEnd = ifelse(contrast %in% comp.order[1:2], 
                          plotting.order,
                          ifelse(contrast %in% comp.order[3:4], 
                                 plotting.order + dodge.width,
                                 plotting.order + dodge.width) ) ) %>%
    ungroup()
  
  # Re-plot with significance bars
  # Ignore Warning message, "Ignoring unknown aesthetics: annotations";
  # this occurs because there is only one value in symbol. That's what I
  # intended!
  r <- q + geom_signif(data = df.thsd,
                       aes(annotations = symbol), manual = T,
                       xmin = df.thsd$xStart, xmax = df.thsd$xEnd, 
                       y_position = df.thsd$y, 
                       annotations = df.thsd$symbol,
                       show.legend = F, color = "black",
                       tip_length = rep(0.04, nrow(df.thsd)),
                       vjust = 0, size = 1, textsize = 0)
  
  # Different plot dimensions are desired for each data subset, because
  # each has a different number of compounds to plot
  plot.width <- length(unique(df$Compound)) * 0.833
  plot.height <- ifelse(df.name == "canonical_monomers_and_dimers", 6,
                        ifelse(df.name == "noncanonical_monomers", 5.25, 
                               ifelse(df.name == "noncanonical_dimers", 7, 6)))
  
  # Export as pdf 
  pdf(paste0("enzymes_", df.name, "_plot.pdf"), 
      width = plot.width, height = plot.height)
  plot(r)
  dev.off()
  
  # Close loop
}

## Create one final plot of the legend only
# First, recreate the last plot but with the legend enabled and formatted
s <- ggplot() +
  geom_point(data = df, show.legend = T,
             aes(x = as.factor(Compound), y = Log.AUC, 
                 fill = as.factor(Sample.type)),
             size = 5, pch = 21, color = "black", 
             position = position_dodge(width = box.width)) +
  # Set the plot dimensions and point colors
  coord_cartesian(ylim = y.lims) +
  # Add colors manually, last two values are for background rectangles
  scale_fill_manual(values = color.values,
                    labels = color.labels ) +
  scale_x_discrete(labels = x.labels) +
  scale_y_continuous(breaks = y.breaks,
                     labels = y.labels) +
  # Make it pretty with some custom settings
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(7, "pt"),
        axis.line = element_line(color = "black", size = 1),
        axis.text = element_text(color= "black", size = 12),
        axis.text.x = element_text(size = 7),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 12),
        legend.position = "bottom",
        legend.key = element_blank())

# extract the legend
leg <- get_legend(s)
l <- as_ggplot(leg)

# export pdf
pdf(paste0("enzymes_legend.pdf"), 
    width = 7, height = 2)
plot(l)
dev.off() 



### Save statistical tests for each dataset

# Convert from numeric back to character convention for Compound and 
# Sample.type in the Compound and contrast columns
log.enzyme.thsd <- log.enzyme.thsd %>%
  mutate( Compound = compound.order[Compound],
          contrast =  apply(str_split(contrast, "-", simplify = T), 1, function(x) {
            paste0( sample.order[as.numeric(x[1])], "-", 
                    sample.order[as.numeric(x[2])] ) 
          } ) )
# Export as .csv files
setwd(wd)
write.csv(log.enzyme.thsd, 
          row.names = F,
          file = paste0("tukeysHSD_enzyme_dataset.csv"))









### Generate plots for A. pisum dataset

## Setup some generic plotting parameters

# Define an order and colors for enzyme.substrate that determines the plotting
# order of samples for each compound. Last two samples/colors are for
# background colors
sample.order <- c("mut.wall", 
                  "sonication.Ap", "lysozyme.Ap")
sample.colors <- c("#440154FF", 
                   "gray60", "gray30")
names(sample.colors) <- c(1:length(sample.order))

# Convert the Compound and Sample.type columns into a numeric convention. 
log.Ap.dat <- log.Ap.dat %>%
  mutate( Compound = as.numeric(factor(Compound, levels = compound.order)),
          Sample.type = as.numeric(factor(Sample.type, levels = sample.order)) ) %>%
  # Arrange the dataset first by compound and then by sample type
  arrange( Compound, Sample.type )

# Extract the Peptide order from the current, sorted data frame
peptide.order <- levels(factor(log.Ap.dat$Peptide , 
                               levels = unique(log.Ap.dat$Peptide)))
# Convert Peptide into a numeric convention based on the peptide order 
log.Ap.dat <- log.Ap.dat %>%
  mutate( Peptide = as.numeric(factor(Peptide, levels = peptide.order)) )


## Add any missing sample data to the dataset  
# Each compound/sample pair should have plot points
log.Ap.dat <- log.Ap.dat %>%
  select( Compound, MS.Alias.name, Peptide, Log.AUC ) %>%
  complete( Compound, MS.Alias.name, 
            fill = list( Log.AUC = 0.30103 ) ) %>%
  group_by( Compound, MS.Alias.name ) %>%
  mutate( Sample.type = case_when( MS.Alias.name %in% c("1A", "2A", "3A") ~ grep("mut.wall", sample.order),
                                   MS.Alias.name == "sonic3" ~ grep("sonication.Ap", sample.order),
                                   MS.Alias.name == "lys3" ~ grep("lysozyme.Ap", sample.order) ),
          Peptide = log.Ap.dat$Peptide[match(Compound, log.Ap.dat$Compound)] ) %>%
  ungroup()


## Split the dataset up further into typical monomers & dimers, 
# atypical monomers, atypical dimers, and trimers, to be plotted separately
df.can <- log.Ap.dat %>%
  filter( Compound %in% match(c(can.mono, can.di, can.tri), compound.order) )
df.ncan.mono <- log.Ap.dat %>%
  filter( Compound %in% match(ncan.mono, compound.order) )
df.ncan.di <- log.Ap.dat %>%
  filter( Compound %in% match(ncan.di, compound.order) )


## Loop through the subsetted data frames
dfs <- list(df.can, df.ncan.mono, df.ncan.di, df.can.tri)
names(dfs) <- c("canonical_monomers_and_dimers", "noncanonical_monomers",
                "noncanonical_dimers", "canonical_trimers")

for (i in 1:length(dfs)) {
  df <- dfs[[i]]
  df.name <- names(dfs)[i]
  
  ## Establish some dataset-specific plotting parameters
  # Create custom y-axis labels for the log scale values
  y.scale <- c( floor(min(log.Ap.dat$Log.AUC)):ceiling(max(log.Ap.dat$Log.AUC)) )
  y.breaks <- c(0: max(y.scale))
  y.labels <- c("ND", "1E-6%", "1E-5%", "1E-4%", "0.001%", 
                "0.01%", "0.1%", "1%", "10%", "100%")
  # Define cartesian limits for y-axis
  y.lims <- c(0, max(y.scale))
  # Define x-axis labels
  x.lims <- c(0, length(unique(df$Compound)))
  x.labels <- compound.order[sort(unique(df$Compound))]
  # Define colors
  color.values <- sample.colors[sort(unique(df$Sample.type))]
  color.labels <- sample.order[sort(unique(df$Sample.type))]
  
  # Establish the width of each boxplot. The actual width of each dodged 
  # boxplot will be this value divided by the number of samples plotted
  # at that point.
  box.width <- 0.75
  
  
  ## Differentiate between different series of compounds that all share the
  ## same stem peptide by using alternating white and gray backgrounds.
  ## This can be accomplished by plotting gray rectangles behind the dot plot
  ## at the appropriate coordinates.
  
  # Define the plotting order based on the Compound order of this dataset
  df <- df %>%
    mutate( plotting.order = as.numeric(factor(Compound , 
                                               levels = unique(Compound))) )
  # Determine coordinates and background colors for each peptide series
  bg.pos <- df %>%
    select( Peptide, plotting.order ) %>%
    group_by( Peptide ) %>%
    summarise( xStart = min(plotting.order), xEnd = max(plotting.order) ) %>%
    # Add a column with the peptide series number for this specific subset
    mutate( series.number = row_number() ) %>%
    # Adjust values so that borders occur between x-values
    mutate( xStart = ifelse(xStart == 1, xStart - 1, xStart - 0.5),
            xEnd = ifelse(xEnd == max(xEnd), Inf, xEnd + 0.5) ) %>%
    # Retrieve y values that have already defined
    mutate ( yStart = -Inf, yEnd = max(y.breaks) ) %>%
    # Specify whether this range is white or gray with "A and "B, respectively
    mutate( Background = ifelse( series.number %% 2 != 0, "A", "B" ) ) %>%
    ungroup()
  
  
  ## Generate plots
  # Calc.M on x-axis, Log.AUC on y-axis, Sample.type is color, and
  # Compound is x-axis label
  q <- ggplot() +
    # First plot the background sections as rectangles
    geom_rect(data = bg.pos, show.legend = F, 
              aes(xmin = xStart, xmax = xEnd, ymin = yStart, ymax = yEnd, 
                  fill = Background)) +
    # Then plot the data points, converting Calc.M and Sample.type to factors 
    geom_point(data = df, show.legend = F,
               aes(x = as.factor(Compound), y = Log.AUC, 
                   fill = as.factor(Sample.type)),
               size = 5, pch = 21, color = "black", 
               position = position_dodge(width = box.width)) +
    # Set the plot dimensions and point colors
    coord_cartesian(ylim = y.lims) +
    # Add colors manually, last two values are for background rectangles
    scale_fill_manual(values = c(color.values, NA, "gray90"),
                      labels = c(color.labels, NA, NA) ) +
    scale_x_discrete(labels = x.labels) +
    scale_y_continuous(breaks = y.breaks,
                       labels = y.labels) +
    # Make it pretty with some custom settings
    annotation_logticks(sides = "l", scaled = T, color = "black") +
    theme(plot.title = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 1, colour = "black"),
          axis.ticks.length = unit(7, "pt"),
          axis.line.x = element_line(color = "black", size = 1),
          axis.text = element_text(color= "black", size = 18),
          axis.text.x = element_text(size = 16, angle = 45, 
                                     hjust = 1, vjust = 1),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 0, 0, 80), "pt") )
  
  # Different plot dimensions are desired for each data subset, because
  # each has a different number of compounds to plot
  plot.width <- length(unique(df$Compound)) * 0.833
  plot.height <- ifelse(df.name == "canonical_monomers_and_dimers", 6,
                        ifelse(df.name == "noncanonical_monomers", 5.25, 
                               ifelse(df.name == "noncanonical_dimers", 6.75, 6)))
  
  # Export as pdf with a plot spacer to the left of the plot so that
  # long x-axis labels don't get lost
  pdf(paste0("Ap_", df.name, "_plot.pdf"), 
      width = plot.width, height = plot.height)
  plot(q)
  dev.off()
}


## Create one final plot of the legend only
# First, recreate the last plot but with the legend enabled and formatted
s <- ggplot() +
  geom_point(data = df, show.legend = T,
             aes(x = as.factor(Compound), y = Log.AUC, 
                 fill = as.factor(Sample.type)),
             size = 5, pch = 21, color = "black", 
             position = position_dodge(width = box.width)) +
  # Set the plot dimensions and point colors
  coord_cartesian(ylim = y.lims) +
  # Add colors manually, last two values are for background rectangles
  scale_fill_manual(values = color.values,
                    labels = color.labels ) +
  scale_x_discrete(labels = x.labels) +
  scale_y_continuous(breaks = y.breaks,
                     labels = y.labels) +
  # Make it pretty with some custom settings
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(7, "pt"),
        axis.line = element_line(color = "black", size = 1),
        axis.text = element_text(color= "black", size = 12),
        axis.text.x = element_text(size = 7),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 12),
        legend.position = "bottom",
        legend.key = element_blank())

# extract the legend
leg <- get_legend(s)
l <- as_ggplot(leg)

# export pdf
pdf(paste0("Ap_legend.pdf"), 
    width = 7, height = 2)
plot(l)
dev.off()

