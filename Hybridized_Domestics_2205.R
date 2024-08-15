library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)
library(radiator)

#### Read in data ####
# Read in 2111 genetic data
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in 2111 metadata
set.seed(27)
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) %>% 
  #filter(Cohort == "Domestic") %>% 
  mutate(Group = case_when(Cohort == "Domestic" ~ sample(1:2, n(), replace = TRUE)))

Samples_2111 %>% count(Group)
  
# Subset by random grouping
Data_2111@pop <- as_factor(Samples_2111$Group)

Group_1 <- popsub(Data_2111, "1")

Group_2 <- popsub(Data_2111, "2")

#### Hybridize domestic fish to create the prototypical "St. Croix Falls domestic" brook trout population ####
hybridize(Group_1, Group_2,
          n = 100,
          pop = "St.Croix_HD",
          res.type = c("genind"),
          hyb.label = "HD") %>% 
  tidy_genind() %>% 
  write_genepop(genepop.header = "St.Croix hybridized domestic fish (n = 100)",
                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/St.Croix_HD")

# Will have to manually fill pop slot when reading in the genepop file as shown below:

St.Croix_HD <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/St.Croix_HD_genepop.gen", 
                            ncode = 3L, 
                            quiet = FALSE)

temp <- tibble(Hybrid_type = rep("St.Croix_HD", 500))
St.Croix_HD@pop <- as_factor(temp$Hybrid_type)
