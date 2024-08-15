library(tidyverse)
library(readxl)
library(miscTools)

#### Create genepop file for converted CE pops ####

# Read in CE population data
CE_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location)

# Read in capillary electrophoresis data and tidy
CE_genotypes_all <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  select(!c("Location", "LocationCode", "SampleCode"))

colnames(CE_genotypes_all) <- CE_genotypes_all %>% 
  colnames() %>% 
  str_replace_all(c("\\." = "_", "-" = "_"))

# Remove loci we aren't interested in AND APPLY CE TO AMPLICON GENOTYPE CONVERSIONS
CE_genotypes <- CE_genotypes_all %>%
  select(-c(sfo115___27, sfo115___28, SfoB52, SfoB52_b, SFO_12, SFO_12_b, SFOC86, SFOC86_b, SFO_18, SFO_18_b)) %>% 
  mutate(L_SFOC113 = case_when(L_SFOC113 != 0 ~ (L_SFOC113 - 79), .default = L_SFOC113),
         L_SFOC113_b = case_when(L_SFOC113_b != 0 ~ (L_SFOC113_b - 79), .default = L_SFOC113_b),
         L_SFOC24 = case_when(L_SFOC24 != 0 ~ (L_SFOC24 - 78), .default = L_SFOC24),
         L_SFOC24_b = case_when(L_SFOC24_b != 0 ~ (L_SFOC24_b - 78), .default = L_SFOC24_b),
         L_SFOC28 = case_when(L_SFOC28 != 0 ~ (L_SFOC28 - 115), .default = L_SFOC28),
         L_SFOC28_b = case_when(L_SFOC28_b != 0 ~ (L_SFOC28_b - 115), .default = L_SFOC28_b),
         L_SFOC88 = case_when(L_SFOC88 != 0 ~ (L_SFOC88 - 114), .default = L_SFOC88),
         L_SFOC88_b = case_when(L_SFOC88_b != 0 ~ (L_SFOC88_b - 114), .default = L_SFOC88_b),
         SfoC38 = case_when(SfoC38 != 0 ~ (SfoC38 - 58), .default = SfoC38),
         SfoC38_b = case_when(SfoC38_b != 0 ~ (SfoC38_b - 58), .default = SfoC38_b),
         SfoD75 = case_when(SfoD75 != 0 ~ (SfoD75 - 115), .default = SfoD75),
         SfoD75_b = case_when(SfoD75_b != 0 ~ (SfoD75_b - 115), .default = SfoD75_b),
         SfoD91 = case_when(SfoD91 != 0 ~ (SfoD91 - 137), .default = SfoD91),
         SfoD91_b = case_when(SfoD91_b != 0 ~ (SfoD91_b - 137), .default = SfoD91_b),
         .keep = "unused")
         
genotype_data_1 <- CE_pop_data %>% 
  arrange(WaterbodyName) %>% # This arrange() is very important to ensure pop delimiters align with pop_counts later on
  select(SampleID) %>% 
  left_join(CE_genotypes, by = "SampleID")

colnames(genotype_data_1) <- genotype_data_1 %>% 
  colnames() %>% 
  str_replace_all(c("\\." = "_", "-" = "_"))

#### Correct names for each locus ####
locus_names <- colnames(genotype_data_1) %>%
  as_tibble() %>%
  slice(-1) %>% 
  rename("Locus" = value) %>% 
  filter(!str_detect(Locus, "_b"))

#### Change MegaSat notation to be missing genepop calls (000) ####
genotype_data_2 <- genotype_data_1 %>%
  select(-SampleID) %>%
  mutate(across(everything(), ~replace(., . ==  0, "000")),
         across(everything(), ~replace(., . ==  "X" , "000")),
         across(everything(), ~replace(., . ==  "Unscored" , "000")),
         across(everything(), ~str_pad(., 3, pad = "0")))

#### Unite the alleles for each locus ####
united_alleles <- genotype_data_1 %>% 
  select(SampleID) %>%   
  mutate(SampleID = paste(genotype_data_1$SampleID, ","))

# 14 because 7 final selected loci
odds <- seq(1, 14, by = 2)

for(i in odds){
  z <- i + 1 
  locus <- colnames(genotype_data_2)[i]
  united_locus <- genotype_data_2 %>% select(i, z) %>% unite({{locus}}, sep = "")
  united_alleles[, ncol(united_alleles) + 1] <- united_locus
}

united_alleles

################### Construct the converted CE Genepop file ################################

genotype_matrix <- as.matrix(united_alleles)

#### Make Genepop header ####
file_date <- format(Sys.time(), "%Y%m%d@%H%M") # date and time
header <- paste("Genepop file format", "MCGL2205", file_date)

#### List of locus names separated by commas ####
locus_names_2 <- paste(locus_names$Locus, collapse = ",")

#### Generate appropriate "Pop" rows ####
# Pop label that will separate each population
pop_line <- c("Pop", rep("", ncol(genotype_matrix)-1))

#### Count the number of individuals in each population ####
pop_counts <- data.frame(Counts = count(CE_pop_data, WaterbodyName))

#### Add a column totalling the cumulative sum ####
pop_counts <- pop_counts %>% mutate(Sum = cumsum(pop_counts$Counts.n))

#### Insert a Pop row between each population ####
for (i in 1:nrow(pop_counts)){
  # i is the row number and increases by 1 after each iteration to compensate
  # for the extra row being inserted each run through the loop
  pop.row <- rep(NA, nrow(pop_counts))
  pop.row[i] <- pop_counts$Sum[i] + i
  genotype_matrix <- insertRow(genotype_matrix, pop.row[i], pop_line)
}

genotype_matrix <- genotype_matrix[-nrow(genotype_matrix), ]

# Insert title, locus and pop rows at the beginning
genotype_matrix <- insertRow(genotype_matrix, 1, c(header, rep("", ncol(genotype_matrix)-1 )))
genotype_matrix <- insertRow(genotype_matrix, 2, c(locus_names_2, rep("", ncol(genotype_matrix)-1 )))
genotype_matrix <- insertRow(genotype_matrix, 3, pop_line)

# Export file
write.table(genotype_matrix, file = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/CE_genepop_converted.gen",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)