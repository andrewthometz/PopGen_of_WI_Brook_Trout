library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)

#### Read in data ####
# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in 2111 genetic data for reference comparisons
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  arrange(Cohort)  %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab)))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

Data_2111@pop <- as_factor(Samples_2111$Cohort)

# Join 2205 and 2111 data
temporary <- popsub(Data_2111, sublist = c("D_2018", "D_2019", "D_2020"))
temporary@pop <- factor(rep("St.Croix domestic", nrow(temporary@tab)))

Joined_2205_2111 <- repool(Data_2205, temporary) 

#### First create a distance matrix to be supplied ####
# hierfstat package (This distance matrix didn't seem to agree with the AMOVA... different size?)
hierfstat_2205 <- genind2hierfstat(Data_2205)
dist_matrix <- as.matrix(genet.dist(hierfstat_2205, method = "Dch"))

# adegenet package (This distance matrix appears to work better)
genpop_list <- as.genpop(Data_2205$tab)
dist_matrix_2 <- dist.genpop(genpop_list)

#### Analysis of Molecular Variance ####
# Fill strata slot of genind object
Data_2205@strata <- Samples_2205

# Link that describes missing data arguments: https://grunwaldlab.github.io/poppr/reference/missingno.html
amova_result <- poppr.amova(Data_2205, 
                            hier = ~ HUC_2 / HUC_4 / HUC_6 / HUC_8 / WaterbodyName,
                            #dist = dist_matrix_2,
                            within = FALSE,
                            missing = "genotype",
                            cutoff = 0.5)
amova_result

amova_signif <- randtest(amova_result) # Test for significance (takes ~ 10 min to run)
plot(amova_signif)

#### Re-run with St.Croix domestics included ####
# Fill strata slot of genind object
strata <- Samples_2111 %>% 
  select(-WaterbodyName) %>% 
  mutate(WaterbodyName = recode(Cohort,
                                "D_2018" = "St.Croix domestic",
                                "D_2019" = "St.Croix domestic",
                                "D_2020" = "St.Croix domestic")) %>% 
  filter(WaterbodyName == "St.Croix domestic")

Joined_2205_2111@strata <- Samples_2205 %>% 
  bind_rows(strata)

amova_2_result <- poppr.amova(Joined_2205_2111, 
                              hier = ~ WaterbodyName,
                              #dist = ,
                              within = FALSE,
                              missing = "genotype",
                              cutoff = 0.5)
amova_2_result

amova_2_signif <- randtest(amova_2_result) # Test for significance (Takes ~ 2 min to run)
plot(amova_2_signif)

