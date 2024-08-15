library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)

#### Read in Master brook trout genepop file ####
UNIFIED_BKT <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/UNIFIED_BKT_genepop.gen",
                            ncode = 3L,
                            quiet = FALSE)

#### Read in metadata ####
# 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls domestic",
         HUC_2 = "Hatchery",
         HUC_4 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_4, HUC_2)

# 2205 metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  #filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  #arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) %>% 
  mutate(Data_source = "Thometz")

# Erdman
Erdman_samples <- read_excel("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx") %>%
  mutate(Data_source = "Erdman")

# Bind the metadata
All_metadata <- Samples_2205 %>% 
  bind_rows(Erdman_samples) %>% 
  select(SampleID, WBIC, WaterbodyName, HUC_4, HUC_8, Data_source, Latitude, Longitude) %>% 
  filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) 

fctLevelOrder <- All_metadata %>% 
  distinct(WaterbodyName) %>% 
  arrange(WaterbodyName)

All_metadata <- All_metadata %>% 
  mutate(WaterbodyName = fct_relevel(WaterbodyName, fctLevelOrder$WaterbodyName))

# Assign pop slot
UNIFIED_BKT@pop <- as_factor(All_metadata$WaterbodyName)

#### Run pairwise Gst Nei ####
gstMatrix <- pairwise_Gst_Nei(UNIFIED_BKT)

gst_tidy <- gstMatrix %>% 
  tidy()

gst_tidy %>% write_csv("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Analyses/Gst_Fst/Gst_4loci_UNIFIED.csv")

heatmap <- gst_tidy %>% 
  mutate(distance = round(distance, digits = 2)) %>% 
  ggplot(aes(x = item1, y = item2, fill = distance)) +
  geom_tile(color = "grey") +
  scale_fill_gradient(low = "blue", 
                      high = "red", 
                      space = "Lab", 
                      name = bquote("Nei's"~G[ST]),
                      breaks = c(round(min(gst_tidy$distance), 2),
                                 round(mean(gst_tidy$distance), 2),
                                 round(max(gst_tidy$distance), 2))) +
  scale_y_discrete(limits = rev(rownames(as.matrix(gstMatrix)))) +
  scale_x_discrete(position = "bottom") +
  labs(x = "",
       y = "") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -60, 
                                   vjust = 0.5, 
                                   size = 6, 
                                   hjust = 0.01,
                                   color = "black"),
        axis.text.y = element_text(size = 6,
                                   color = "black"))

ggsave(filename = "GST_heatmap_4loci_UNIFIED.pdf",
       plot = heatmap,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 12,
       width = 14,
       units = "in")

ggsave(filename = "GST_heatmap_4loci_UNIFIED.png",
       plot = heatmap,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 12,
       width = 14,
       units = "in")
