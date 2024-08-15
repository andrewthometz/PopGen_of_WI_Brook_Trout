library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(radiator)
library(RColorBrewer)
library(ggrepel)

library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

#### Read in Amplicon genepop file ####
Amp_gen <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Amplicon_genepop.gen", 
                        ncode = 3L, 
                        quiet = FALSE)

# Read in CE vs Amp population data
Amp_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  filter(SampleID %in% rownames(Amp_gen@tab)) %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location) %>% 
  arrange(match(SampleID, rownames(Amp_gen@tab)))

# Fill pop slot
Amp_gen@pop <- as_factor(Amp_pop_data$WaterbodyName)

#### Read in data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  arrange(Cohort) %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) %>% 
  mutate(WaterbodyName = case_when(str_detect(WaterbodyName, "St. Croix") ~ "Hatchery",
                                   .default = WaterbodyName))

# Fill 2111 pop slot and subset to just domestic fish
Data_2111@pop <- as_factor(Samples_2111$WaterbodyName)

Domestics <- popsub(Data_2111, 
                    sublist = "Hatchery",
                    drop = FALSE)

# Filter to CE vs Amp loci
Domestics <- Domestics[loc = locNames(Amp_gen)]

# Read in 2205 genetic data for reference comparisons
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# # Filter to CE vs Amp loci
# temp_gen <- Data_2205[loc = locNames(Amp_gen)]
# 
# # Write and read it back in to work around repooling error
# temp_gen %>% 
#   tidy_genind() %>% 
#   write_genepop(genepop.header = "59 Survey pops with just 7 CE vs Amp loci",
#                 filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/2205_pops_7_loci")

#Data_2205_7_loci <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/2205_pops_7_loci.gen", 
#                                 ncode = 3L, 
#                                 quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205_7_loci@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205_7_loci@tab))) %>% 
  select(WaterbodyName, SampleID, WBIC, HUC_8, HUC_6, HUC_4, HUC_2)

Data_2205_7_loci@pop <- as_factor(Samples_2205$WaterbodyName)

# Repool
All_dat_amplicon <- repool(Data_2205_7_loci, Domestics, Amp_gen)

# Create df with all SampleID's and Waterbodies
centroid_df_amp <- Samples_2111 %>% 
  bind_rows(Samples_2111) %>% 
  bind_rows(Amp_pop_data) %>% 
  select(SampleID, WaterbodyName)

#### Run PCA with amplicon data ####
pca_amp_data <- tab(All_dat_amplicon, freq = TRUE, NA.method = "mean")

pca_amp_result <- dudi.pca(pca_amp_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

pop_names_1 <- All_dat_amplicon@pop %>% 
  as_tibble() %>% 
  rename(pop = value)

pca_df_1 <- pca_amp_result$li %>% 
  as_tibble() %>% 
  mutate(Population = pop_names_1$pop)

# Plot
centroids_pop_1 <- pca_df_1 %>% 
  select(Population) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ Population, pca_df_1, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_amp_plot <- pca_df_1 %>% 
  ggplot(aes(x = Axis1, y = Axis2)) + 
  #geom_point(alpha = 0.4) + 
  #stat_ellipse() + 
  geom_point(data = centroids_pop_1, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_pop_1, 
                  aes(label = Population), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

ggsave(filename = "PCA_amplicon.pdf",
       plot = pca_amp_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 12,
       width = 12,
       units = "in")

######################## Do this again with converted CE data ##################################

#### Read in converted CE genepop file ####
CE_converted <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/CE_genepop_converted.gen", 
                             ncode = 3L, 
                             quiet = FALSE)

# Read in CE population data
CE_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  filter(SampleID %in% rownames(CE_converted@tab)) %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location) %>% 
  arrange(match(SampleID, rownames(CE_converted@tab)))

# Fill 2111 pop slot and subset to just domestic fish
CE_converted@pop <- as_factor(CE_pop_data$WaterbodyName)

# Repool
All_dat_converted <- repool(Data_2205_7_loci, Domestics, CE_converted)

# Create df with all SampleID's and Waterbodies
centroid_df_CE <- Samples_2111 %>% 
  bind_rows(Samples_2111) %>% 
  bind_rows(CE_pop_data) %>% 
  select(SampleID, WaterbodyName)

#### Run PCA with amplicon data ####
pca_CE_data <- tab(All_dat_converted, freq = TRUE, NA.method = "mean")

pca_CE_result <- dudi.pca(pca_CE_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

pop_names_2 <- All_dat_converted@pop %>% 
  as_tibble() %>% 
  rename(pop = value)

pca_df_2 <- pca_CE_result$li %>% 
  as_tibble() %>% 
  mutate(Population = pop_names_2$pop)

# Plot
centroids_pop_2 <- pca_df_2 %>% 
  select(Population) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ Population, pca_df_2, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_CE_plot <- pca_df_2 %>% 
  ggplot(aes(x = Axis1, y = Axis2)) + 
  #geom_point(alpha = 0.4) + 
  #stat_ellipse() + 
  geom_point(data = centroids_pop_2, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_pop_2, 
                  aes(label = Population), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

ggsave(filename = "PCA_converted_CE.pdf",
       plot = pca_CE_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 12,
       width = 12,
       units = "in")