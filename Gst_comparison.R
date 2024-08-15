library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)
library(xlsx)

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

Amp_pop_data$WaterbodyName <- paste(Amp_pop_data$WaterbodyName, "_Amp", sep = "")

# Fill pop slot
Amp_gen@pop <- as_factor(Amp_pop_data$WaterbodyName)

# Read in converted CE data
CE_converted <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/CE_genepop_converted.gen", 
                             ncode = 3L, 
                             quiet = FALSE)

# Read in CE population data
CE_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  filter(SampleID %in% rownames(CE_converted@tab)) %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location) %>% 
  arrange(match(SampleID, rownames(CE_converted@tab)))

CE_pop_data$WaterbodyName <- paste(CE_pop_data$WaterbodyName, "_CE", sep = "")

# Fill 2111 pop slot and subset to just domestic fish
CE_converted@pop <- as_factor(CE_pop_data$WaterbodyName)

# Repool
CE_and_amp <- repool(CE_converted, Amp_gen)

################ run it and plot it

gstMatrix <- pairwise_Gst_Nei(CE_and_amp)

#write.xlsx(as.matrix(gstMatrix), file = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/CE_Amp_Gst_matrix.xlsx") 

heatmap <- gstMatrix %>% 
  tidy() %>% 
  ggplot(aes(x = item1, y = item2, fill = distance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", 
                      high = "red", 
                      space = "Lab", 
                      name = bquote(G[ST])) +
  scale_y_discrete(limits = rev(rownames(as.matrix(gstMatrix)))) +
  scale_x_discrete(position = "bottom") +
  xlab("") +
  ylab("") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -60, 
                                   vjust = 0.5, 
                                   size = 18, 
                                   hjust = 0.05,
                                   color = "black"),
        axis.text.y = element_text(size = 16,
                                   color = "black"),
        legend.title = element_text(size = 36), # change legend title font size
        legend.text = element_text(size = 28),
        legend.key.size = unit(0.5, 'in'))

ggsave(filename = "GST_heatmap_CE_Amp.pdf",
       plot = heatmap,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 10,
       width = 14,
       units = "in")

#### Check for missing data in amp genotypes ####
Amp_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location)

# Read in capillary electrophoresis data and tidy
Amp_genotypes_all <- read_delim("X:/2201_BKT_msat_conversion/Sfon_2201-001_paired_uSat_output/Sfon_2201-001_Genotype.txt") %>% 
  rename(SampleID = Sample_idx1_idx2)

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

test <- genotype_data_1 %>% 
  left_join(Amp_pop_data)

#### Make scatter plot comparing CE vs Amp comparisons ####
gstMatrix_Amp <- pairwise_Gst_Nei(Amp_gen) %>% 
  tidy() %>% 
  rename(Amp_distance = distance)

gstMatrix_CE <- pairwise_Gst_Nei(CE_converted) %>% 
  tidy() %>% 
  rename(CE_distance = distance)

CE_Amp_tib <- tibble(gstMatrix_Amp$Amp_distance, gstMatrix_CE$CE_distance) %>% 
  mutate(CE_distance = round(gstMatrix_CE$CE_distance, 3),
         Amp_distance = round(gstMatrix_Amp$Amp_distance, 3),
         .keep = "used") %>% 
  mutate(difference = Amp_distance - CE_distance)

scatter <- CE_Amp_tib %>% 
  ggplot(aes(x = CE_distance, y = Amp_distance)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  xlim(0, 0.1) +
  ylim(0, 0.1) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "CE distance",
       y = "Amp distance") +
  theme_classic()

ggsave(filename = "Scatterplot_CE_Amp.pdf",
       plot = scatter,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 10,
       width = 14,
       units = "in")

