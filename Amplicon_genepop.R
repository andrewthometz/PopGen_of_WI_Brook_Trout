library(tidyverse)
library(readxl)
library(miscTools)

#### Create genepop file for non-converted CE pops ####

# Read in CE population data
Amp_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location) 

# Read in capillary electrophoresis data and tidy
Amp_genotypes_all <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_Genotype.txt") %>% 
  rename(SampleID = Sample_idx1_idx2) %>% 
  mutate(SampleID = case_when(SampleID == "17-03009" ~ "17-03051",
                              SampleID == "17-03010" ~ "17-03052",
                              SampleID == "17-03011" ~ "17-03053",
                              SampleID == "17-03012" ~ "17-03054",
                              SampleID == "17-03013" ~ "17-03055",
                              SampleID == "17-03014" ~ "17-03056",
                              SampleID == "17-03015" ~ "17-03057",
                              SampleID == "17-03016" ~ "17-03058",
                              
                              SampleID == "17-03051" ~ "17-03009",
                              SampleID == "17-03052" ~ "17-03010",
                              SampleID == "17-03053" ~ "17-03011",
                              SampleID == "17-03054" ~ "17-03012",
                              SampleID == "17-03055" ~ "17-03013",
                              SampleID == "17-03056" ~ "17-03014",
                              SampleID == "17-03057" ~ "17-03015",
                              SampleID == "17-03058" ~ "17-03016",
                              .default = SampleID)) # this mutate corrects swapping of columns 8 and 9 in lab plate

colnames(Amp_genotypes_all) <- Amp_genotypes_all %>% 
  colnames() %>% 
  str_replace_all(c("\\." = "_", "-" = "_"))

# Remove loci we aren't interested in
Amp_genotypes <- Amp_genotypes_all %>%
  select(SampleID,
         L_SFOC113,
         L_SFOC113_b,
         L_SFOC24,
         L_SFOC24_b,
         L_SFOC28,
         L_SFOC28_b,
         L_SFOC88,
         L_SFOC88_b,
         SfoC38,
         SfoC38_b,
         SfoD75,
         SfoD75_b,
         SfoD91,
         SfoD91_b)

genotype_data_1 <- Amp_pop_data %>% 
  arrange(WaterbodyName) %>% # This arrange() is very important to ensure pop delimiters align with pop_counts later on
  select(SampleID) %>% 
  left_join(Amp_genotypes, by = "SampleID")

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

################### Construct the non-converted CE Genepop file ################################

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
pop_counts <- data.frame(Counts = count(Amp_pop_data, WaterbodyName))

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
write.table(genotype_matrix, file = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Amplicon_genepop.gen",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)