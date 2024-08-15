library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(poppr)
library(radiator)

#### Create new genepop file to incorporate the 4 pops from project 2113 into my original 59 pops ####

# Original 59 pops
Data_2205_temp <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_genepop_59pops.gen", 
                               ncode = 3L, 
                               quiet = FALSE)



# Read in 2113 genepop file for additional 4 pops
Data_2113_temp <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2113_genepop_4pops.gen", 
                               ncode = 3L, 
                               quiet = FALSE)

# Subset loci down to the 68 loci from my project
Data_2113 <- Data_2113_temp[loc = locNames(Data_2205_temp)]

# Repool the first two genepops
Data_2205 <- repool(Data_2113, Data_2205_temp)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

# Write master genepop file
Data_2205 %>% 
  tidy_genind() %>% 
  write_genepop(genepop.header = "59 pops from 2205 + 4 pops from 2113 @ 68 loci",
                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_All_63pops")
