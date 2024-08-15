library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(poppr)
library(radiator)

#### Create new genepop file to incorporate random subsample of 30 St. Croix domestic fish for genetic structure analyses ####

#### Read in data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) %>% 
  mutate(WaterbodyName = case_when(str_detect(WaterbodyName, "St. Croix") ~ "St. Croix Falls domestic",
                                   .default = WaterbodyName))

# Fill 2111 pop slot and subset to just domestic fish
Data_2111@pop <- as_factor(Samples_2111$WaterbodyName)

Domestics <- popsub(Data_2111, 
                    sublist = "St. Croix Falls domestic",
                    drop = FALSE)

# subsample 30 fish
Domestics_sample <- Domestics[sample(x = 450, size = 30)]

# Read in 2205 genetic data for reference comparisons
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_All_63pops.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

# Repool
Joined_2205_2111 <- repool(Data_2205, Domestics_sample)

# Write master genepop file
Joined_2205_2111 %>% 
  tidy_genind() %>% 
  write_genepop(genepop.header = "All 63 pops + random subsample of 30 St. Croix Falls domestic fish @ 68 loci",
                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics")
