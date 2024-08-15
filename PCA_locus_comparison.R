library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(radiator)
library(ggrepel)

#### Read in Master brook trout genepop file ####
UNIFIED_BKT <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/UNIFIED_BKT.gen",
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

Erdman_samples %>% 
  filter(WaterbodyName %in% Samples_2205$WaterbodyName) %>% 
  count(WaterbodyName)

# Bind the metadata
All_metadata <- Samples_2205 %>% 
  bind_rows(Erdman_samples) %>% 
  select(SampleID, WBIC, WaterbodyName, HUC_4, HUC_8, Data_source, Latitude, Longitude) %>% 
  filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) 

test <- All_metadata %>% count(WaterbodyName)

# Assign pop slot
UNIFIED_BKT@pop <- as_factor(All_metadata$WaterbodyName)

# Create df with all SampleID's and Waterbodies
centroid_df <- All_metadata %>% 
  select(SampleID, WaterbodyName)

#### Run PCA at each locus ####
# L_SFOC113
L_SFOC113_gen <- UNIFIED_BKT[loc = "L_SFOC113"]

L_SFOC113_data <- tab(L_SFOC113_gen, freq = TRUE, NA.method = "mean")

pca_113_result <- dudi.pca(L_SFOC113_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

summary(pca_113_result)

pca_df_113 <- pca_113_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName,
         HUC_4 = All_metadata$HUC_4,
         HUC_8 = All_metadata$HUC_8)

pop_points_113 <- pca_113_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName) %>%
  group_by(WaterbodyName) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

huc_points_113 <- pca_113_result$li %>% 
  as_tibble() %>% 
  mutate(HUC_4 = All_metadata$HUC_4) %>%
  group_by(HUC_4) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

pca_113_plot <- pca_df_113 %>% 
  ggplot(aes(x = Axis1, y = Axis2)) + 
  #stat_ellipse(geom = "polygon", 
  #             alpha = 0.25,
  #             aes(color = HUC_4)) + 
  geom_point(alpha = 0.25) +
  geom_point(data = pop_points_113, 
             size = 1, 
             alpha = 0.5) +
  geom_point(data = huc_points_113, 
             size = 3, 
             alpha = 0.75) + 
  geom_text_repel(data = huc_points_113, 
                  aes(label = HUC_4), 
                  fontface = "bold", 
                  size = 3, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  labs(x = "Axis 1 (39.6%)",
       y = "Axis 2 (17.1%)",
       title = "L_SFOC113") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

ggsave(filename = "PCA_L_SFOC113.pdf",
       plot = pca_113_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 12,
       units = "in")

# L_SFOC28
L_SFOC28_gen <- UNIFIED_BKT[loc = "L_SFOC28"]

L_SFOC28_data <- tab(L_SFOC28_gen, freq = TRUE, NA.method = "mean")

pca_28_result <- dudi.pca(L_SFOC28_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

summary(pca_28_result)

pca_df_28 <- pca_28_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName,
         HUC_4 = All_metadata$HUC_4,
         HUC_8 = All_metadata$HUC_8)

pop_points_28 <- pca_28_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName) %>%
  group_by(WaterbodyName) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

huc_points_28 <- pca_28_result$li %>% 
  as_tibble() %>% 
  mutate(HUC_4 = All_metadata$HUC_4) %>%
  group_by(HUC_4) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

pca_28_plot <- pca_df_28 %>% 
  ggplot(aes(x = Axis1, y = Axis2)) + 
  #stat_ellipse(geom = "polygon", 
  #             alpha = 0.25,
  #             aes(color = HUC_4)) + 
  geom_point(alpha = 0.25) +
  geom_point(data = pop_points_28, 
             size = 1, 
             alpha = 0.5) +
  geom_point(data = huc_points_28, 
             size = 3, 
             alpha = 0.75) + 
  geom_text_repel(data = huc_points_28, 
                  aes(label = HUC_4), 
                  fontface = "bold", 
                  size = 3, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  labs(x = "Axis 1 (34.6%)",
       y = "Axis 2 (27.3%)",
       title = "L_SFOC28") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

ggsave(filename = "PCA_L_SFOC28.pdf",
       plot = pca_28_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 12,
       units = "in")

# L_SFOC88
L_SFOC88_gen <- UNIFIED_BKT[loc = "L_SFOC88"]

L_SFOC88_data <- tab(L_SFOC88_gen, freq = TRUE, NA.method = "mean")

pca_88_result <- dudi.pca(L_SFOC88_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

summary(pca_88_result)

pca_df_88 <- pca_88_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName,
         HUC_4 = All_metadata$HUC_4,
         HUC_8 = All_metadata$HUC_8)

pop_points_88 <- pca_88_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName) %>%
  group_by(WaterbodyName) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

huc_points_88 <- pca_88_result$li %>% 
  as_tibble() %>% 
  mutate(HUC_4 = All_metadata$HUC_4) %>%
  group_by(HUC_4) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

pca_88_plot <- pca_df_88 %>% 
  ggplot(aes(x = Axis1, y = Axis2)) + 
  #stat_ellipse(geom = "polygon", 
  #             alpha = 0.25,
  #             aes(color = HUC_4)) + 
  geom_point(alpha = 0.25) +
  geom_point(data = pop_points_88, 
             size = 1, 
             alpha = 0.5) +
  geom_point(data = huc_points_88, 
             size = 3, 
             alpha = 0.75) + 
  geom_text_repel(data = huc_points_88, 
                  aes(label = HUC_4), 
                  fontface = "bold", 
                  size = 3, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  labs(x = "Axis 1 (53.6%)",
       y = "Axis 2 (23.2%)",
       title = "L_SFOC88") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

ggsave(filename = "PCA_L_SFOC88.pdf",
       plot = pca_88_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 12,
       units = "in")

# Sfoc38
Sfoc38_gen <- UNIFIED_BKT[loc = "SfoC38"]

Sfoc38_data <- tab(Sfoc38_gen, freq = TRUE, NA.method = "mean")

pca_38_result <- dudi.pca(Sfoc38_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

summary(pca_38_result)

pca_df_38 <- pca_38_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName,
         HUC_4 = All_metadata$HUC_4,
         HUC_8 = All_metadata$HUC_8)

pop_points_38 <- pca_38_result$li %>% 
  as_tibble() %>% 
  mutate(WaterbodyName = All_metadata$WaterbodyName) %>%
  group_by(WaterbodyName) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

huc_points_38 <- pca_38_result$li %>% 
  as_tibble() %>% 
  mutate(HUC_4 = All_metadata$HUC_4) %>%
  group_by(HUC_4) %>% 
  summarize(Axis1 = mean(Axis1),
            Axis2 = mean(Axis2))

pca_38_plot <- pca_df_38 %>% 
  ggplot(aes(x = Axis1, y = Axis2)) + 
  #stat_ellipse(geom = "polygon", 
  #             alpha = 0.25,
  #             aes(color = HUC_4)) + 
  geom_point(alpha = 0.25) +
  geom_point(data = pop_points_38, 
             size = 1, 
             alpha = 0.5) +
  geom_point(data = huc_points_38, 
             size = 3, 
             alpha = 0.75) + 
  geom_text_repel(data = huc_points_38, 
                  aes(label = HUC_4), 
                  fontface = "bold", 
                  size = 3, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  labs(x = "Axis 1 (59.3%)",
       y = "Axis 2 (25.0%)",
       title = "SfoC38") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

ggsave(filename = "PCA_SfoC38.pdf",
       plot = pca_38_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 12,
       units = "in")
