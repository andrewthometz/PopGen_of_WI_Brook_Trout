library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(poppr)
library(radiator)

# Read in my CE converted genepop for Erdman's data, prep for repooling
Erdman_gen <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_converted.gen",
                           ncode = 3L,
                           quiet = FALSE)

# Read in gen file including all 63 pops plus 30 random St. Croix Falls domestics
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls domestic",
         HUC_2 = "Hatchery",
         HUC_4 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_4, HUC_2)

Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

Data_2205 <- Data_2205[loc = locNames(Erdman_gen)]

# Have to write and read it back in to work around repooling error (also for 4 loci vs 68 loci comparisons)
#Data_2205 %>% 
#  tidy_genind() %>% 
#  write_genepop(genepop.header = "63 survey pops + 30 St. Croix Falls domestics @ 4 loci",
#                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/2205_63pops_4loci")

Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/2205_63pops_4loci.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in my CE converted genepop for Erdman's data, prep for repooling
Erdman_gen <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_converted.gen",
                           ncode = 3L,
                           quiet = FALSE)

Erdman_samples <- read_excel("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx") %>% 
  select(SampleID, Pop, WaterbodyName, Latitude, Longitude) %>% 
  filter(SampleID %in% rownames(Erdman_gen@tab)) %>% 
  arrange(match(SampleID, rownames(Erdman_gen@tab)))

Erdman_gen@pop <- as_factor(Erdman_samples$WaterbodyName)

Erdman_gen <- popsub(Erdman_gen, exclude = popNames(Data_2205))

# Repool the data and write out master genepop file (4 loci)
repool(Erdman_gen, Data_2205) %>% 
  tidy_genind() %>% 
  write_genepop(genepop.header = "Erdman pops + 63 pops from 2205 + 30 St. Croix domestics @ 4 loci",
                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/UNIFIED_BKT")
