library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(ggh4x)
library(patchwork)

library(sf)
library(ggrepel)
library(ggspatial)

# DAPC tests a hypothesis, PCA does not

# DAPC guidelines from Thia 2022:
# n.da = k groups (should be determined a priori, # of sample pops)
# n.pca must be =< k-1 (only k-1 PCs are biologically informative)

#### Read in data ####
# Read in 2205 genetic data
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

DAPC_plot1 <- DAPC_1_probs_longer %>% 
  filter(HUC_2 == "Upper Mississippi Region") %>% 
  ggplot(aes(x = SampleID, y = Probability, fill = Cluster)) +
  geom_col(show.legend = FALSE) + 
  facet_nested(cols = vars(HUC_8, WaterbodyName),
               switch = "x",
               nest_line = element_line(linewidth = 1, lineend = "round"),
               solo_line = TRUE,
               resect = unit(0.05, "in"),
               scales = "free", 
               space = "free") +
  labs(x = "", 
       y = "Admixture\nproportion",
       title = "Upper Mississippi Region (HUC 2)") +
  scale_y_continuous(expand = c(0, 0),
                     position = "left",
                     breaks = seq(0, 1, by = 0.5)) +
  scale_fill_manual(values = K7_colors) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0))

DAPC_plot2 <- DAPC_1_probs_longer %>% 
  filter(HUC_2 == "Great Lakes Region" |
         HUC_2 == "Hatchery") %>% 
  ggplot(aes(x = SampleID, y = Probability, fill = Cluster)) +
  geom_col(show.legend = TRUE) + 
  facet_nested(cols = vars(HUC_8, WaterbodyName),
               switch = "x",
               nest_line = element_line(linewidth = 1, lineend = "round"),
               solo_line = TRUE,
               resect = unit(0.05, "in"),
               scales = "free", 
               space = "free") +
  labs(x = "", 
       y = "Admixture\nproportion",
       title = "Great Lakes Region (HUC 2)",
       fill = "Cluster (DAPC)") +
  scale_y_continuous(expand = c(0, 0),
                     position = "left",
                     breaks = seq(0, 1, by = 0.5)) +
  scale_fill_manual(values = K7_colors) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1))

DAPC_plots <- DAPC_plot1 / DAPC_plot2

ggsave(filename = "DAPC_str_plot.pdf",
       plot = DAPC_plots,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/DAPC",
       height = 10,
       width = 14,
       units = "in")

ggsave(filename = "DAPC_str_plot.png",
       plot = DAPC_plots,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/DAPC",
       height = 10,
       width = 14,
       units = "in")

############################
#### Plot them on a map ####

# Read in necessary shape files
HUC8_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp")

HUC2_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Major_Basins/Major_Basins.shp")

WMU_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Water_Management_Units/Water_Management_Units.shp")

# Prep admixture df and lat long df for mapmixture function
DAPC_mapmixture <- DAPC_1_probs %>% 
  left_join(Samples_2205) %>% 
  select(WaterbodyName, SampleID, K1, K2, K3, K4, K5, K6, K7) %>% 
  filter(WaterbodyName != "St. Croix Falls Strain")

# These are intentionally incorrect, revised for ease of viewing
Lats_Longs <- Samples_2205 %>% 
  select(WaterbodyName, Latitude, Longitude) %>% 
  filter(WaterbodyName != "St. Croix Falls Strain") %>% 
  distinct() %>% 
  mutate(Longitude = case_when(WaterbodyName == "Swan Creek" ~ -91.1,
                               WaterbodyName == "Marshall Creek - West Branch" ~ -90.7,
                               WaterbodyName == "Unnamed trib to Dell Creek (b)" ~ -89.97,
                               WaterbodyName == "Fourmile Creek" ~ -89.4,
                               WaterbodyName == "Lunch Creek" ~ -89.55,
                               WaterbodyName == "Lowery Creek" ~ -90.1,
                               WaterbodyName == "Knapp Creek" ~ -90.75,
                               WaterbodyName == "Tagatz Creek" ~ -89.75,
                               WaterbodyName == "Plover River" ~ -89.4,
                               WaterbodyName == "Marshall Creek" ~ -90.45,
                               WaterbodyName == "Flume Creek" ~ -89.15, .default = Longitude),
         Latitude = case_when(WaterbodyName == "Bruce Creek" ~ 43.99,
                              WaterbodyName == "Marshall Creek" ~ 43.34,
                              WaterbodyName == "Little Willow Creek" ~ 45.7,
                              WaterbodyName == "Alvin Creek" ~ 45.91, .default = Latitude))

# Plot using mapmixture function
DAPC_map <- mapmixture(admixture_df = DAPC_mapmixture, 
                       coords_df = Lats_Longs,
                       basemap = HUC2_shp,
                       boundary = c(xmin = -93.5, 
                                    xmax = -86.5, 
                                    ymin = 42, 
                                    ymax = 47.5),
                       cluster_names = c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7"),
                       cluster_cols = K7_colors,
                       pie_size = 0.25,
                       pie_border = 0.05,
                       land_colour = "grey80",
                       sea_colour = NA,
                       arrow_size = 2,
                       arrow_position = "tr",
                       scalebar_size = 1,
                       scalebar_position = "tr",
                       axis_title_size = 8,
                       axis_text_size = 6,
                       plot_title = "DAPC (K = 7)",
                       plot_title_size = 10) +
    theme_classic() +
    guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(fill = "Cluster") +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.text = element_text(size = 6))

ggsave(filename = "DAPC_map.pdf",
       plot = DAPC_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/DAPC",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "DAPC_map.png",
       plot = DAPC_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/DAPC",
       height = 5,
       width = 5,
       units = "in")

#############################################
#### Try running DAPC as ordination plot ####
ind_clusters <- DAPC_1_probs_longer %>% 
  group_by(SampleID) %>% 
  filter(Probability == max(Probability))

ind_coords <- DAPC_1_optimal$ind.coord %>% 
  data.frame() %>% 
  rownames_to_column(var = "SampleID")

centroid_coords <- DAPC_1_optimal$grp.coord %>% 
  data.frame() %>% 
  rownames_to_column(var = "Cluster") %>% 
  mutate(Cluster = case_when(Cluster == 4 ~ "K1",
                             Cluster == 3 ~ "K2",
                             Cluster == 2 ~ "K3",
                             Cluster == 1 ~ "K4",
                             Cluster == 6 ~ "K5",
                             Cluster == 5 ~ "K6",
                             Cluster == 7 ~ "K7"),
         Cluster = str_replace(Cluster, "K", ""),
         Cluster = str_replace(Cluster, "1", "1 (St. Croix Falls Strain)"),
         Cluster = fct_relevel(Cluster, c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7")))
                               
# df 1 and 2
ord_plot_1 <- ind_coords %>% 
  ggplot(aes(x = LD1, y = LD2, color = ind_clusters$Cluster)) +
  geom_point(alpha = 0.25) +
  stat_ellipse(alpha = 0.75,
               level = 0.95,
               linewidth = 0.75) +
  geom_point(data = centroid_coords,
             aes(color = Cluster),
             size = 5) +
  geom_label_repel(data = centroid_coords,
            aes(label = Cluster,
                color = Cluster),
            size = 4,
            #force = 1.5,
            #force_pull = 2,
            label.padding = 0.15
            ) +
  labs(x = "Discriminant function 1",
       y = "Discriminant function 2",
       title = "(A)   DAPC discriminant functions 1 & 2"
       #color = "Genetic cluster\n(DAPC)"
       ) +
  scale_color_manual(values = K7_colors) +
  theme_classic() +
  theme(legend.position = "none")
  
# df 2 and 3
ord_plot_2 <- ind_coords %>% 
  ggplot(aes(x = LD2, y = LD3, color = ind_clusters$Cluster)) +
  geom_point(alpha = 0.25) +
  stat_ellipse(alpha = 0.75,
               level = 0.95,
               linewidth = 0.75) +
  geom_point(data = centroid_coords,
             aes(color = Cluster),
             size = 5) +
  geom_label_repel(data = centroid_coords,
            aes(label = Cluster,
                color = Cluster),
            size = 4,
            #force = 1.5,
            #force_pull = 2,
            label.padding = 0.15
            ) +
  labs(x = "Discriminant function 2",
       y = "Discriminant function 3",
       title = "(B)   DAPC discriminant functions 2 & 3"
       #color = "Genetic cluster\n(DAPC)"
       ) +
  scale_color_manual(values = K7_colors) +
  theme_classic() +
  theme(legend.position = "none")

ord_plot <- ord_plot_1 / ord_plot_2

ggsave(filename = "DAPC_ordplot.pdf",
       plot = ord_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/DAPC",
       height = 9,
       width = 9,
       units = "in")

ggsave(filename = "DAPC_ordplot.png",
       plot = ord_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/DAPC",
       height = 9,
       width = 9,
       units = "in")

## Create tree to visualize cluster relatedness (and ensure correct color assignment for plots) ##
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
    values = K7_colors) +
  xlim(0, 0.2) + # This can help make tree fit
  labs(title = "DAPC clusters (K = 7)")

ggsave(filename = "DAPC_Cluster_Tree.pdf",
       plot = tree_1,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 4,
       width = 7,
       units = "in")

ggsave(filename = "DAPC_Cluster_Tree.png",
       plot = tree_1,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 4,
       width = 7,
       units = "in")


