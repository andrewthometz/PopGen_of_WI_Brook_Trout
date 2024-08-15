library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(RColorBrewer)
library(patchwork)
library(ggh4x)


library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

# Create a 3-panel figure including all three cluster trees

################ DAPC K = 7 ###########################################

Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Prep 2111 data to work with plotting
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain",
         HUC_2 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_2)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

#### DAPC as a means of understanding genetic structure ####
# Run a DAPC to determine the number of clusters in the data #  https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
set.seed(27)
clusters <- find.clusters.genind(Data_2205, 
                                 max.n.clust = 20,
                                 n.pca = 300,
                                 n.clust = 7
) 

# Select 300 to retain all PCs. Takes a bit to run
# Selecting 7 clusters, as that's where the lowest BIC scores begin to level out

table(Data_2205$pop, clusters$grp)

DAPC_1 <- dapc.genind(Data_2205, 
                      pop = clusters$grp, # 7 clusters
                      n.pca = nPop(Data_2205) - 1,
                      n.da = nPop(Data_2205)
)

summary(DAPC_1) # Low assignment proportions indicate admixture, high indicate clear-cut clusters

# Use A-score to find optimal number of PCs
A_score_1 <- optim.a.score(DAPC_1)   # Suggests 7 PCs is optimal

# Re-run with optimal number of PCs
DAPC_1_optimal <- dapc.genind(Data_2205, 
                              pop = clusters$grp, # 7 clusters
                              n.pca = 7, 
                              n.da = nPop(Data_2205))
summary(DAPC_1_optimal)

##########################
#### Plot the results ####

# Create custom color palette
brewer.pal(n = 7, name = "Set1")
K7_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3", "#A65628")

# Plot cluster membership probabilities
DAPC_1_probs <- DAPC_1_optimal$posterior %>%
  as_tibble(rownames = "SampleID") %>% 
  rename(K1 = `4`,
         K2 = `3`,
         K3 = `2`,
         K4 = `1`,
         K5 = `6`,
         K6 = `5`,
         K7 = `7`)

DAPC_1_probs_longer <- DAPC_1_probs %>%
  pivot_longer(cols = 2:8, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205) %>% 
  mutate(Cluster = str_replace(Cluster, "K", ""),
         Cluster = str_replace(Cluster, "1", "1 (St. Croix Falls Strain)"),
         Cluster = fct_relevel(Cluster, c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7")),
         HUC_8 = fct_relevel(HUC_8, "Hatchery", after = Inf))

# Filter to fish with at least 75% assignment to a given cluster
DAPC_probs_filtered <- DAPC_1_probs_longer %>% 
  filter(Probability >= 0.75) %>% 
  select(SampleID, Cluster)

# Revise Samples_2205 to make popsub possible
Samples_revised <- Samples_2205 %>% 
  left_join(DAPC_probs_filtered) %>% 
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
tree_DAPC <- ggtree(Tree_data, 
                 aes(color = Cluster), 
                 size = 1,
                 show.legend = FALSE) +
  geom_tiplab(show.legend = FALSE) +
  #geom_treescale(color = "black",
  #               linesize = 1) +
  geom_text(aes(label = Bootstraps),
            hjust = -0.25,
            size = 2,
            show.legend = FALSE) +
  scale_colour_manual(#name = "Cluster", # Can use this or scale_color_discrete()
    na.value = "black",
    values = K7_colors) +
  xlim(0, 0.2) + # This can help make tree fit
  labs(title = "(A)   DAPC clusters (K = 7)")

############################## K = 6 #############################

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
K6 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Final_run_2205/AssProbs_CleanedUp/K6_AssProbs_CleanedUp.txt") %>% 
  select(-c(n, percent_miss)) %>% 
  mutate(C6 = as.numeric(C6)) %>% 
  rename(K1 = C4,
         K2 = C6,
         K3 = C3,
         K4 = C1,
         K5 = C5,
         K6 = C2)

K6_longer <- K6 %>% 
  pivot_longer(cols = 2:7, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205) %>% 
  mutate(Cluster = str_replace(Cluster, "K", ""),
         Cluster = str_replace(Cluster, "1", "1 (St. Croix Falls Strain)"),
         Cluster = fct_relevel(Cluster, c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6")),
         HUC_8 = fct_relevel(HUC_8, "Hatchery", after = Inf))

# Filter to fish with at least 75% assignment to a given cluster
K6_filtered <- K6_longer %>% 
  filter(Probability >= 0.75) %>% 
  select(SampleID, Cluster)

# Revise Samples_2205 to make popsub possible
Samples_revised <- Samples_2205 %>% 
  left_join(K6_filtered) %>% 
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
brewer.pal(n = 6, name = "Set1")
K6_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3")

tree_K6 <- ggtree(Tree_data, 
                 aes(color = Cluster), 
                 size = 1,
                 show.legend = FALSE) +
  geom_tiplab(show.legend = FALSE) +
  #geom_treescale(color = "black",
  #               linesize = 1) +
  geom_text(aes(label = Bootstraps),
            hjust = -0.25,
            size = 2,
            show.legend = FALSE) +
  scale_colour_manual(#name = "Cluster", # Can use this or scale_color_discrete()
    na.value = "black",
    values = K6_colors) +
  xlim(0, 0.2) + # This can help make tree fit
  labs(title = "(B)   STRUCTURE clusters (K = 6)")

############################## K = 9 ##########################################

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

tree_K9 <- ggtree(Tree_data, 
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
  labs(title = "(C)   STRUCTURE clusters (K = 9)")

################################# Join them together and save #########################

combined_tree <- tree_DAPC / tree_K6 / tree_K9

ggsave(filename = "Combined_Cluster_Tree_3panel.pdf",
       plot = combined_tree,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 8,
       width = 7,
       units = "in")

ggsave(filename = "Combined_Cluster_Tree_3panel.png",
       plot = combined_tree,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 8,
       width = 7,
       units = "in")

