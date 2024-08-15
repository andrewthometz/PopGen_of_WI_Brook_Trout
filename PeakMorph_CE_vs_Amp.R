library(tidyverse)
library(ggforce)

#### Replicate of the PeakMorph_2111 script but for samples where CE and Amp calls are in disagreement ####

### Function to make length_dist files tidy ###

length_dist_tidy <- function(x){
  
  BKT_ID <- deparse(substitute(x))
  
  output <- x %>% 
    mutate(uSat_locus = str_replace_all(Microsatellite, "-", "_"), .keep = "unused") %>% 
    select(-sum) %>%
    pivot_longer(-c(uSat_locus, scores), names_to = "Length", values_to = "Read_count") %>% 
    add_column(SampleID = BKT_ID) %>% 
    mutate(SampleID = str_remove_all(SampleID, "BKT_")) %>% 
    mutate(Length = as.numeric(Length))
  
  print(output)
}

### Read in 12 BKT ### These are the 12 fish with the greatest number of CE vs Amp disagreements

BKT_17_03053 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03053.txt")
BKT_17_03009 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03009.txt")
BKT_17_03011 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03011.txt")

BKT_17_03051 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03051.txt")
BKT_17_03052 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03052.txt")
BKT_17_03010 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03010.txt")

BKT_17_03015 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03015.txt")
BKT_17_03016 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03016.txt")
BKT_17_03058 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03058.txt")
BKT_17_03057 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03057.txt")

BKT_17_03012 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03012.txt")
BKT_17_03054 <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_length_distribution/Genotype_17-03054.txt")

### Run each fish's length_dist file through my tidying function while joining them ###

All_LD_data <- bind_rows(length_dist_tidy(BKT_17_03053),
                         length_dist_tidy(BKT_17_03009),
                         length_dist_tidy(BKT_17_03011),
                         length_dist_tidy(BKT_17_03051),
                         length_dist_tidy(BKT_17_03052),
                         length_dist_tidy(BKT_17_03010),
                         length_dist_tidy(BKT_17_03015),
                         length_dist_tidy(BKT_17_03016),
                         length_dist_tidy(BKT_17_03058),
                         length_dist_tidy(BKT_17_03057),
                         length_dist_tidy(BKT_17_03012),
                         length_dist_tidy(BKT_17_03054)) %>% 
  filter(uSat_locus == "L_SFOC113" |
         uSat_locus == "L_SFOC24" |
         uSat_locus == "L_SFOC28" |
         uSat_locus == "L_SFOC88" |
         uSat_locus == "SFOC86" |
         uSat_locus == "SFO_18" |
         uSat_locus == "SfoC38" |
         uSat_locus == "SfoD75" |
         uSat_locus == "SfoD91")

All_LD_data %>% 
  summarize(n_loci = n_distinct(uSat_locus))

All_LD_data %>% 
  #filter(Read_count > 0) %>% 
  group_by(uSat_locus, SampleID) %>% 
  #slice_max(Read_count, n = 10) %>% 
  filter(Read_count >= max(Read_count)*0.25) %>% 
  ungroup() %>% 
  summarize(x = n_distinct(uSat_locus))

### Make peak morphology plots, one locus and 12 BKT per page ### Saves directly as pdf

pdf("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/PeakMorph_CE_vs_Amp.pdf", paper = "a4r", width = 11, height = 9)

ProgressBar <- txtProgressBar(min = 0, max = 9, style = 3)

for(i in 1:9){
  
print(All_LD_data %>% 
        group_by(uSat_locus, SampleID) %>% 
        #slice_max(Read_count, n = 8) %>% ### This seems to determine run time (4 ~ 40min, 10 ~ 1.5hrs)
        filter(Read_count >= max(Read_count)*0.05) %>% 
        ggplot(aes(x = Length, y = Read_count)) +
        geom_col(fill = "seagreen") +
        geom_text(aes(label = Length), 
                  hjust = 0.5, 
                  vjust = 1.5, 
                  colour = "black") +
        ylab("Read count") +
        xlab("Allele length") +
        scale_y_continuous(expand = c(0,0)) +
        #scale_x_continuous(breaks = c(0:200), expand = c(0,0)) +
        theme_classic() +
        facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(SampleID), as.factor(scores)), 
                            nrow = 3,
                            ncol = 4,
                            page = i, 
                            scales = "free"))
  Sys.sleep(0.1)
  setTxtProgressBar(ProgressBar, i)
}
close(ProgressBar)
dev.off()

### Single page of peak morphs for testing plot changes ###

All_LD_data %>% 
  #filter(Read_count > 0) %>% 
  group_by(uSat_locus, SampleID) %>% 
  #slice_max(Read_count, n = 10) %>% 
  filter(Read_count >= max(Read_count)*0.05) %>% 
  ggplot(aes(x = Length, y = Read_count)) +
  geom_col(fill = "seagreen") +
  geom_text(aes(label = Length), hjust = 0.5, vjust = 1.5, colour = "black") +
  ylab("Read count") +
  xlab("Allele length") +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(breaks = c(0:200), expand = c(0,0)) +
  theme_classic() +
  facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(SampleID), as.factor(scores)), 
                      nrow = 3,
                      ncol = 4,
                      page = 49, 
                      scales = "free")
