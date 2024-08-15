library(tidyverse)
library(readxl)
library(poppr)
library(adegenet)
library(ggrepel)

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

# Read in 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) %>% 
  mutate(WaterbodyName = case_when(str_detect(WaterbodyName, "St. Croix") ~ "St. Croix Falls domestic",
                                   .default = WaterbodyName))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

Data_2111@pop <- as_factor(Samples_2111$WaterbodyName)

# Join 2205 and 2111 data
domestics <- popsub(Data_2111, "St. Croix Falls domestic")

Joined_2205_2111 <- repool(Data_2205, domestics)

#### Run a PCA that just shows centroid (see mcglPCA function code) and color centroids by huc level ####
pca_data <- tab(Data_2205, freq = TRUE, NA.method = "mean")

pca_result <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

pop_names <- Data_2205@pop %>% 
  as_tibble() %>% 
  rename(pop = value)

pca_df <- pca_result$li %>% 
  as_tibble() %>% 
  mutate(Population = Samples_2205$WaterbodyName,
         HUC_2 = Samples_2205$HUC_2,
         HUC_4 = Samples_2205$HUC_4,
         HUC_6 = Samples_2205$HUC_6,
         HUC_8 = Samples_2205$HUC_8)

# HUC 2
centroids_HUC2 <- pca_df %>% 
  select(HUC_2) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ HUC_2, pca_df, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = HUC_2)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_HUC2, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_HUC2, 
                  aes(label = HUC_2), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

# HUC 2 axis 2 and 3
centroids_HUC2_23 <- pca_df %>% 
  select(HUC_2) %>% 
  right_join(aggregate(cbind(Axis2, Axis3) ~ HUC_2, pca_df, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df %>% 
  ggplot(aes(x = Axis2, y = Axis3, color = HUC_2)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_HUC2_23, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_HUC2_23, 
                  aes(label = HUC_2), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 3") + 
  xlab("Axis 2") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

# HUC 4
centroids_HUC4 <- pca_df %>% 
  select(HUC_4) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ HUC_4, pca_df, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = HUC_4)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_HUC4, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_HUC4, 
                  aes(label = HUC_4), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

# HUC 6
centroids_HUC6 <- pca_df %>% 
  select(HUC_6) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ HUC_6, pca_df, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = HUC_6)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_HUC6, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_HUC6, 
                  aes(label = HUC_6), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

# HUC 8
centroids_HUC8 <- pca_df %>% 
  select(HUC_8) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ HUC_8, pca_df, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = HUC_8)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_HUC8, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_HUC8, 
                  aes(label = HUC_8), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

#### PCA of populations with no stocking record ####

# Filter to find HUC12s with no bkt stocking history (17 pops)
# Read in WDNR stocking database #
Stocking_data <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Stocking_histories_2205.csv")

never_stocked <- Samples_2205 %>% 
  filter(!(HUC_12 %in% Stocking_data$HUC_12))

# Make vector of pop names that haven't been stocked
unstocked_list <- never_stocked %>% 
  distinct(WaterbodyName)

Unstocked_genind <- popsub(Data_2205, sublist = unstocked_list$WaterbodyName)

####### This plot will go in my publication, needs to be cleaned up #######
####### Justifies pooling the native pops together when hybridizing #######
# Run PCA
pca_unstocked <- tab(Unstocked_genind, freq = TRUE, NA.method = "mean")

pca_unstocked_result <- dudi.pca(pca_unstocked, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

pop_names_2 <- Unstocked_genind@pop %>% 
  as_tibble() %>% 
  rename(pop = value)

pca_df_2 <- pca_unstocked_result$li %>% 
  as_tibble() %>% 
  mutate(Population = never_stocked$WaterbodyName,
         HUC_2 = never_stocked$HUC_2,
         HUC_4 = never_stocked$HUC_4,
         HUC_6 = never_stocked$HUC_6,
         HUC_8 = never_stocked$HUC_8)

# HUC 2
centroids_HUC2 <- pca_df_2 %>% 
  select(HUC_2) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ HUC_2, pca_df_2, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df_2 %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = HUC_2)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_HUC2, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_HUC2, 
                  aes(label = HUC_2), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

# Population level
centroids_pop <- pca_df_2 %>% 
  select(Population) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ Population, pca_df_2, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df_2 %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = Population)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_pop, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_pop, 
                  aes(label = Population), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

# Population level including St.Croix hybrids
St.Croix_HD <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/St.Croix_HD_genepop.gen", 
                            ncode = 3L, 
                            quiet = FALSE)

temp <- tibble(WaterbodyName = rep("St.Croix_HD", 500))

St.Croix_HD@pop <- as_factor(temp$WaterbodyName)

Joined_genind <- repool(Unstocked_genind, St.Croix_HD)

fish_df <- never_stocked %>% 
  bind_rows(temp)

# Run pca again with pooled data
pca_unstocked <- tab(Joined_genind, freq = TRUE, NA.method = "mean")

pca_unstocked_result <- dudi.pca(pca_unstocked, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

pop_names_2 <- Joined_genind@pop %>% 
  as_tibble() %>% 
  rename(pop = value)

pca_df_2 <- pca_unstocked_result$li %>% 
  as_tibble() %>% 
  mutate(Population = fish_df$WaterbodyName,
         HUC_2 = fish_df$HUC_2,
         HUC_4 = fish_df$HUC_4,
         HUC_6 = fish_df$HUC_6,
         HUC_8 = fish_df$HUC_8)

centroids_pop <- pca_df_2 %>% 
  select(Population) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ Population, pca_df_2, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df_2 %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = Population)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_pop, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_pop, 
                  aes(label = Population), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

#### PCA comparing each hybrid group #### 
St.Croix_HD <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/St.Croix_HD_genepop.gen", 
                            ncode = 3L, 
                            quiet = FALSE)

temp <- tibble(Hybrid_type = rep("St.Croix_HD", 500))
St.Croix_HD@pop <- as_factor(temp$Hybrid_type)

Hybrid_natives <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/Hybridized_natives_genepop.gen",
                               ncode = 3L, 
                               quiet = FALSE)

temp <- tibble(Hybrid_type = rep("Hybrid_native", 500))
Hybrid_natives@pop <- as_factor(temp$Hybrid_type)

Hybrids_only <- repool(Hybrid_natives, St.Croix_HD)

# Run pca again with pooled hybrid data
pca_hybrids <- tab(Hybrids_only, freq = TRUE, NA.method = "mean")

pca_hybrids_result <- dudi.pca(Hybrids_only, center = TRUE, scale = FALSE, nf = 4, scannf = FALSE)

Pop_names <- Hybrids_only@pop %>% 
  as_tibble() %>% 
  rename(pop = value)

pca_df_hybrids <- pca_hybrids_result$li %>% 
  as_tibble() %>% 
  mutate(Hybrid_type = Pop_names$pop)

centroids_type <- pca_df_hybrids %>% 
  select(Hybrid_type) %>% 
  right_join(aggregate(cbind(Axis1, Axis2) ~ Hybrid_type, pca_df_hybrids, mean)) %>% 
  distinct(.keep_all = TRUE)

pca_df_hybrids %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = Hybrid_type)) + 
  geom_point(alpha = 0.4) + 
  stat_ellipse() + 
  geom_point(data = centroids_type, 
             size = 6, 
             alpha = 0.8) + 
  geom_text_repel(data = centroids_type, 
                  aes(label = Hybrid_type), 
                  fontface = "bold", 
                  size = 5, force = 5, 
                  force_pull = 0.1, 
                  max.overlaps = 100) + 
  ylab("Axis 2") + 
  xlab("Axis 1") + 
  theme_classic(base_size = 18) + 
  theme(legend.position = "none")

#####


