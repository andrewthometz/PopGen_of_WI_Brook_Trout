library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(RColorBrewer)

library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

# Prep 2111 data to work with plotting
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain",
         HUC_2 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_2)

# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

# K = 9 STRUCTURE run
K9 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Final_run_2205/AssProbs_CleanedUp/K9_AssProbs_CleanedUp.txt") %>% 
  select(-c(n, percent_miss)) %>% 
  mutate(C9 = as.numeric(C9)) %>% 
  rename(K1 = C2,
         K2 = C4,
         K3 = C3,
         K4 = C9,
         K5 = C6,
         K6 = C5,
         K7 = C7,
         K8 = C8,
         K9 = C1)

K9_longer <- K9 %>% 
  pivot_longer(cols = 2:10, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205) %>% 
  mutate(Cluster = str_replace(Cluster, "K", ""),
         Cluster = str_replace(Cluster, "1", "1 (St. Croix Falls Strain)"),
         Cluster = fct_relevel(Cluster, c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7", "8", "9")),
         HUC_8 = fct_relevel(HUC_8, "Hatchery", after = Inf))

# Filter to fish with at least 75% assignment to a given cluster
K9_filtered <- K9_longer %>% 
  filter(Probability >= 0.75) %>% 
  select(SampleID, Cluster)

# Revise Samples_2205 to make popsub possible
Samples_revised <- Samples_2205 %>% 
  left_join(K9_filtered) %>% 
  mutate(Cluster = case_when(is.na(Cluster) ~ "sub_75", .default = Cluster))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_revised$Cluster)

Data_2205_filtered <- popsub(Data_2205, exclude = "sub_75")

# Build initial tree (creates phylo object)
Phylo_tree <- aboot(Data_2205_filtered,
                    strata = Data_2205_filtered@pop,
                    distance = "nei.dist",
                    cutoff = 1,
                    tree = "nj") # Do "nj" instead of default "upgma" to make dendrogram

# Turn phylo object into tibble to add huc data
tree_tibble <- Phylo_tree %>% 
  as_tibble() %>% 
  mutate(Cluster = case_when(label %in% Samples_revised$Cluster ~ label, .default = NA),
         Bootstraps = case_when(!(label %in% Samples_revised$Cluster) ~ label, .default = NA))

# Convert tibble into treedata object for ggtree plotting
Tree_data <- as.treedata(tree_tibble)

# Plot treedata object using ggtree (dendrogram)
brewer.pal(n = 9, name = "Set1")
K9_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3", "#A65628", "#F781BF", "#999999")

tree_1 <- ggtree(Tree_data, 
                 aes(color = Cluster), 
                 size = 1,
                 show.legend = FALSE) +
  geom_tiplab(show.legend = FALSE) +
  geom_treescale(color = "black",
                 linesize = 1) +
  geom_text(aes(label = Bootstraps),
            hjust = -0.25,
            size = 2,
            show.legend = FALSE) +
  scale_colour_manual(#name = "Cluster", # Can use this or scale_color_discrete()
                      na.value = "black",
                      values = K9_colors) +
  xlim(0, 0.2) + # This can help make tree fit
  labs(title = "STRUCTURE clusters (K = 9)")

ggsave(filename = "K9_Cluster_Tree.pdf",
       plot = tree_1,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 4,
       width = 7,
       units = "in")

ggsave(filename = "K9_Cluster_Tree.png",
       plot = tree_1,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 4,
       width = 7,
       units = "in")

#######################################################
#### Build full tree and color code by cluster now ####
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

# Build initial tree (creates phylo object)
Phylo_tree_2 <- aboot(Data_2205,
                      strata = Data_2205@pop,
                      distance = "nei.dist",
                      cutoff = 50,
                      tree = "nj")

temp_info <- K9_longer %>% 
  select(WaterbodyName, Cluster) %>% 
  distinct()

pop_clust <- K9_longer %>% 
  group_by(WaterbodyName, Cluster) %>% 
  summarize(mean_prob = mean(Probability)) %>% 
  group_by(WaterbodyName) %>% 
  filter(mean_prob == max(mean_prob)) %>% 
  left_join(temp_info) %>%
  #filter(WaterbodyName != "St. Croix Falls Strain") %>% 
  rename(label = WaterbodyName)

# Turn phylo object into tibble to add huc data
tree_tibble_2 <- Phylo_tree_2 %>% 
  as_tibble() %>% 
  left_join(pop_clust)

my_breaks <- c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7", "8", "9")

# Convert tibble into treedata object for ggtree plotting
Tree_data_2 <- as.treedata(tree_tibble_2)

# Plot treedata object using ggtree (dendrogram)
#K9_colors_corrected <- c("#4DAF4A", "#377EB8","#FF7F00", "yellow3",  "#E41A1C", "#984EA3", "#999999", "#F781BF", "#A65628")

brewer.pal(n = 9, name = "Set1")
K9_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3", "#A65628", "#F781BF", "#999999")

tree_2 <- ggtree(Tree_data_2, 
                 aes(color = Cluster), 
                 size = 1) +
  geom_tiplab(show.legend = FALSE) +
  geom_treescale(color = "black",
                 linesize = 1) +
  #geom_text(aes(label = Bootstraps),
  #          hjust = -0.25,
  #          size = 3,
  #          show.legend = FALSE) +
  scale_colour_manual(name = "STRUCTURE\nCluster\n(K = 9)",
                      na.value = "black",
                      values = K9_colors,
                      breaks = my_breaks
                      ) +
  xlim(0, 0.25) + # This can help make tree fit
  theme_tree(legend.position = c(0.83, 0.4),
             legend.key.size = unit(0.25, 'in'),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Full_Tree_K9.pdf",
       plot = tree_2,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 10,
       width = 9,
       units = "in")

ggsave(filename = "Full_Tree_K9.png",
       plot = tree_2,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 10,
       width = 9,
       units = "in")

