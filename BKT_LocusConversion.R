library(readxl)
library(tidyverse)
library(adegenet)
library(poppr)

#### Create conversion factor for CE and Amplicon genotype data ####

# Read in final locus selections
Locus_data <- read_excel("X:/2111_F1F2D_BKT/BKT_Locus_Evaluation.xlsx") %>% 
  select(1:9)

hwe_cutoff <- 0.15
naf_cutoff <- 0.2

Final_loci <- Locus_data %>% 
  filter(naf_2205 < naf_cutoff &
           HWE_prop_2205 < hwe_cutoff |
           Locus == "SfoC38" | 
           Locus == "SFOC86",
           Locus != "Salv_1_Di_30_1186") %>% 
  select(Locus)

# Read in capillary electrophoresis data and tidy
CE_genotypes <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  select(!c("Location", "LocationCode", "SampleCode"))

colnames(CE_genotypes) <- CE_genotypes %>% 
  colnames() %>% 
  str_replace_all(c("\\." = "_", "-" = "_"))

# Contains all legacy loci included in final selected loci
CE_tidy <- CE_genotypes %>% 
  pivot_longer(cols = 2:25,
               names_to = "Locus",
               values_to = "CE_allele") %>% 
  filter(Locus != "sfo115___27", 
         Locus != "sfo115___28",
         Locus != "SfoB52",
         Locus != "SfoB52_b",
         Locus != "SFO_12",
         Locus != "SFO_12_b") #%>% 
  #mutate(Locus = str_replace_all(Locus, "_b", ""), .keep = "unused")

# Read in amplicon data and tidy
Amplicon_genotypes <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_Genotype.txt") %>% 
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
                              .default = SampleID))

colnames(Amplicon_genotypes) <- Amplicon_genotypes %>% 
  colnames() %>% 
  str_replace_all(c("\\." = "_", "-" = "_"))

# Contains all legacy loci included in final selected loci
Amp_tidy <- Amplicon_genotypes %>% 
  pivot_longer(cols = 2:183,
               names_to = "Locus",
               values_to = "Amp_allele") %>% 
  #mutate(Locus = str_replace_all(Locus, "_b", ""), .keep = "unused") %>% 
  filter(Locus %in% CE_tidy$Locus,
         SampleID %in% CE_tidy$SampleID) %>% 
  mutate(Amp_allele = as.numeric(Amp_allele))

# Bring in amplicon read depth data
Amplicon_depths <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_coverage.txt") %>% 
  select(-sum) %>% 
  pivot_longer(cols = 2:97,
               names_to = "SampleID",
               values_to = "Depth") %>% 
  mutate(Locus = str_replace_all(Microsatellite, c("\\." = "_", "-" = "_")), .keep = "unused")

# Join CE + Amp genotypes and calculate allele differences
CE_vs_Amp <- CE_tidy %>% 
  left_join(Amp_tidy) %>% 
  drop_na() %>% 
  filter(Amp_allele != 0,
         CE_allele != 0) %>% 
  mutate(Difference = CE_allele - Amp_allele) %>% 
  mutate(Locus = str_replace_all(Locus, "_b", ""), .keep = "unused") %>% 
  left_join(Amplicon_depths)

# Plot
Allele_diffs_count <- CE_vs_Amp %>% 
  filter(Depth > 50) %>% 
  group_by(Locus) %>% 
  count(Difference) %>%
  ggplot(aes(x = Difference, y = n)) +
  geom_col(fill = "seagreen") +
  geom_text(aes(label = Difference)) +
  labs(x = "Allele difference", y = "Number of occurences") +
  facet_wrap(~Locus,
             scales = "free") +
  theme_classic()

ggsave(filename = "Allele_diffs_count.pdf",
       plot = Allele_diffs_count,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 10,
       width = 14,
       units = "in")  

Allele_diffs_prop <- CE_vs_Amp %>% 
  filter(Depth > 50) %>% 
  group_by(Locus) %>% 
  count(Difference) %>% 
  mutate(Proportion = round(n/sum(n), 2)) %>% 
  ggplot(aes(x = Difference, y = Proportion)) +
  geom_col(fill = "royalblue") +
  geom_text(aes(label = Difference)) +
  labs(x = "Allele difference", y = "Proportion of occurences") +
  facet_wrap(~Locus,
             scales = "free") +
  theme_classic()

ggsave(filename = "Allele_diffs_prop.pdf",
       plot = Allele_diffs_prop,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 10,
       width = 14,
       units = "in")  

#### Isolate some samples to run through peak morph script ####
fish_errors <- CE_vs_Amp %>% 
  mutate(error = case_when(Locus == "L_SFOC113" & Difference != 79 ~ "yes",
                           Locus == "L_SFOC113" & Difference == 79 ~ "no",
                           Locus == "L_SFOC24" & Difference != 78 ~ "yes",
                           Locus == "L_SFOC24" & Difference == 78 ~ "no",
                           Locus == "L_SFOC28" & Difference != 115 ~ "yes",
                           Locus == "L_SFOC28" & Difference == 115 ~ "no",
                           Locus == "L_SFOC88" & Difference != 114 ~ "yes",
                           Locus == "L_SFOC88" & Difference == 114 ~ "no",
                           Locus == "SFOC86" & Difference != 40 ~ "yes",
                           Locus == "SFOC86" & Difference == 40 ~ "no",
                           Locus == "SFO_18" & Difference != 39 ~ "yes",
                           Locus == "SFO_18" & Difference == 39 ~ "no",
                           Locus == "SfoC38" & Difference != 58 ~ "yes",
                           Locus == "SfoC38" & Difference == 58 ~ "no",
                           Locus == "SfoD75" & Difference != 115 ~ "yes",
                           Locus == "SfoD75" & Difference == 115 ~ "no",
                           Locus == "SfoD91" & Difference != 137 ~ "yes",
                           Locus == "SfoD91" & Difference == 137 ~ "no")) %>%
  filter(error == "yes") %>% 
  group_by(Locus) %>% 
  count(Locus) %>% 
  rename(n_locus_errors = n) %>% 
  arrange(desc(n_locus_errors)) # Sending these top 12 fish through my peak morph script

#write_csv(fish_errors, file = "X:/2111_F1F2D_BKT/Fish_w_loc_diffs.csv")

#### Conversion error plots ####
pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  select(Location, SampleID) %>% 
  rename(WaterbodyName = Location)

all_pops <- CE_vs_Amp %>% 
  mutate(error = case_when(Locus == "L_SFOC113" & Difference != 79 ~ "yes",
                           Locus == "L_SFOC113" & Difference == 79 ~ "no",
                           Locus == "L_SFOC24" & Difference != 78 ~ "yes",
                           Locus == "L_SFOC24" & Difference == 78 ~ "no",
                           Locus == "L_SFOC28" & Difference != 115 ~ "yes",
                           Locus == "L_SFOC28" & Difference == 115 ~ "no",
                           Locus == "L_SFOC88" & Difference != 114 ~ "yes",
                           Locus == "L_SFOC88" & Difference == 114 ~ "no",
                           Locus == "SFOC86" & Difference != 40 ~ "yes",
                           Locus == "SFOC86" & Difference == 40 ~ "no",
                           Locus == "SFO_18" & Difference != 39 ~ "yes",
                           Locus == "SFO_18" & Difference == 39 ~ "no",
                           Locus == "SfoC38" & Difference != 58 ~ "yes",
                           Locus == "SfoC38" & Difference == 58 ~ "no",
                           Locus == "SfoD75" & Difference != 115 ~ "yes",
                           Locus == "SfoD75" & Difference == 115 ~ "no",
                           Locus == "SfoD91" & Difference != 137 ~ "yes",
                           Locus == "SfoD91" & Difference == 137 ~ "no")) %>% 
  filter(Locus != "SFOC86",
         Locus != "SFO_18") %>% # not using these loci
  left_join(pop_data)

all_pops_plot <- all_pops %>% 
  group_by(SampleID) %>% 
  count(error) %>% 
  left_join(pop_data) %>% 
  ggplot(aes(x = SampleID, y = n, fill = error)) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_col() +
  labs(x = "Sample ID",
       y = "Count") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  facet_wrap(facets = "WaterbodyName",
             scales = "free")

ggsave(filename = "All_fish_errors.pdf",
       plot = all_pops_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 10,
       width = 14,
       units = "in")  

temp_col <- all_pops %>% 
  group_by(Locus) %>% 
  select(WaterbodyName, Locus)

all_locus_plot <- all_pops %>% 
  group_by(Locus, WaterbodyName) %>% 
  count(error) %>% 
  #left_join(pop_data) %>% 
  #bind_cols(temp_col$WaterbodyName) +
  ggplot(aes(x = Locus, y = n, fill = error)) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_col() +
  labs(x = "Locus",
       y = "Count") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  facet_wrap(facets = "WaterbodyName",
             scales = "free")

ggsave(filename = "All_locus_errors.pdf",
       plot = all_locus_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 10,
       width = 14,
       units = "in")  
