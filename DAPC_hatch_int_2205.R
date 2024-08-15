library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(ggh4x)
library(patchwork)

# DAPC tests a hypothesis, PCA does not

# DAPC guidelines from Thia 2022:
# n.da = k groups (should be determined a priori, # of sample pops)
# n.pca must be =< k-1 (only k-1 PCs are biologically informative)

#### Read in data ####
# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_All_63pops.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

# Create HUC df for later
HUCs <- Samples_2205 %>% 
  select(WBIC, WaterbodyName, HUC_2, HUC_4, HUC_6, HUC_8) %>% 
  distinct()

# Read in St.Croix hybrid data
St.Croix_HD <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/St.Croix_HD_genepop.gen", 
                            ncode = 3L, 
                            quiet = FALSE)

temp <- tibble(WaterbodyName = rep("St.Croix_HD", 100))
St.Croix_HD@pop <- as_factor(temp$WaterbodyName)

# Read in native hybrid data
Hybrid_natives <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/Hybridized_natives.gen",
                               ncode = 3L, 
                               quiet = FALSE)

temp <- tibble(WaterbodyName = rep("Hybrid_native", 100))
Hybrid_natives@pop <- as_factor(temp$WaterbodyName)

# Join hybrid data (have to pool all together then popsub to make predict.dapc() function work)
All_fish <- repool(Data_2205, Hybrid_natives, St.Croix_HD)

Hybrid_fish <- popsub(All_fish, c("Hybrid_native", "St.Croix_HD"), drop = FALSE) # drop=F is imperative to retain equal number of variables

Wild_fish <- popsub(All_fish, exclude = c("Hybrid_native", "St.Croix_HD"), drop = FALSE)

#####################################################################################
#### Run DAPC using supplementary individuals to quantify hatchery introgression ####
clusters_hybrids <- find.clusters.genind(Hybrid_fish, 
                                         #max.n.clust = 5,
                                         n.pca = 300,
                                         n.clust = 2)
# Select 300 to retain all PCs
# Selecting 2 clusters, lowest BIC

table(Hybrid_fish$pop, clusters_hybrids$grp)

DAPC_hybrids <- dapc.genind(Hybrid_fish, 
                            #pop = clusters_hybrids$grp,
                            n.clust = 2,
                            n.pca = nPop(Hybrid_fish) - 1,
                            n.da = nPop(Hybrid_fish))

summary(DAPC_hybrids)

#A_score_hybrids <- optim.a.score(DAPC_hybrids) # Not necessary because only 1 PC

# Bring 2205 wild populations in as supplementary individuals to assess hatchery vs native assignments
prediction <- predict.dapc(DAPC_hybrids, newdata = Wild_fish)

prediction$assign
prediction$posterior
prediction$ind.scores

Assignment_probs <- round(prediction$posterior, 2) %>%
  as_tibble(rownames = "SampleID") %>% 
  mutate(WaterbodyName = Data_2205@pop) %>% 
  rename(Hatchery = "St.Croix_HD", Native = "Hybrid_native") %>% 
  pivot_longer(cols = 2:3, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(HUCs) %>% 
  mutate(Cluster = case_when(Cluster == "Hatchery" ~ "St. Croix Falls Strain",
                             .default = Cluster))

#############################################################
#### Quantify St.Croix influence ####
# Overall mean St.Croix influence (Mean assignment probability to hatchery group)
Assignment_probs %>% 
  filter(Cluster == "St. Croix Falls Strain") %>% 
  summarize(overall_mean = round(mean(Probability), 3)*100,
            overall_sd = round(sd(Probability), 3)*100)

# Means by NJ tree groupings
Tree_groupings <- read_csv("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/Trees/Tree_groups.csv")

Group_hatch_ID <- Assignment_probs %>% 
  filter(Cluster == "St. Croix Falls Strain") %>% 
  left_join(Tree_groupings) %>% 
  group_by(Tree_group) %>% 
  summarize(Ave_hatchery_ID = round(mean(Probability), 3),
            SD = round(sd(Probability), 3))

# Mean St.Croix influence (Mean assignment probability to hatchery group)
Ave_hatchery_ID <- Assignment_probs %>% 
  filter(Cluster == "St. Croix Falls Strain") %>% 
  group_by(WaterbodyName) %>% 
  summarize(Ave_hatchery_ID = round(mean(Probability), 3))

# Join together and write as csv
Samples_2205 %>% 
  select(WBIC, WaterbodyName) %>% 
  distinct() %>% 
  left_join(Ave_hatchery_ID) %>% 
  write_csv("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Hatchery_introgression/Hatchery_ID_2205.csv")

# Plot all assignment probabilities
HI_plot_1 <- Assignment_probs %>% 
  filter(HUC_2 == "Upper Mississippi Region") %>% 
  mutate(Cluster = fct_relevel(Cluster, c("Native", "St. Croix Falls Strain"))) %>%
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
  scale_fill_manual(values = c("#1E88E5", "#D81B60")) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0))

HI_plot_2 <- Assignment_probs %>% 
  filter(HUC_2 == "Great Lakes Region" #|
         #HUC_2 == "Hatchery"
         ) %>% 
  mutate(Cluster = fct_relevel(Cluster, c("Native", "St. Croix Falls Strain"))) %>%
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
       fill = "Reference group") +
  scale_y_continuous(expand = c(0, 0),
                     position = "left",
                     breaks = seq(0, 1, by = 0.5)) +
  scale_fill_manual(values = c("#1E88E5", "#D81B60")) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1))

HI_plots <- HI_plot_1 / HI_plot_2
  
ggsave(filename = "Hatch_int_str_plot.pdf",
       plot = HI_plots,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hatch_int",
       height = 10,
       width = 14,
       units = "in")

ggsave(filename = "Hatch_int_str_plot.png",
       plot = HI_plots,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hatch_int",
       height = 10,
       width = 14,
       units = "in")

##################
#### Plotting ####

# scatter(DAPC_hybrids)

# Create density plot for original DAPC model (like scatter(DAPC_hybrids))
LD1_coords <- DAPC_hybrids$ind.coord %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(LD1 = round(LD1, 2),
         Group = Hybrid_fish@pop,
         Group = fct_relevel(Group, c("St.Croix_HD", "Hybrid_native"))) 

density_plot <- LD1_coords %>% 
  ggplot(aes(x = LD1, fill = Group)) +
  geom_density(alpha = 0.6) +
  labs(x = "Discriminant function 1",
       y = "Density",
       title = "(a)",
       fill = "Reference group") +
  scale_fill_manual(values = c("#D81B60", "#1E88E5"),
                    labels = c("St. Croix Falls Strain", "Native")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "top")

# Create figure for reassignment accuracy
hybrid_reassignment_plot <- summary(DAPC_hybrids)$assign.per.pop %>% 
  as_tibble(rownames = "Group") %>% 
  mutate(reassignment = round(value*100, 0)) %>%
  ggplot(aes(x = Group, y = reassignment)) +
  geom_col() +
  coord_flip() +
  labs(x = "Reference group",
       y = "% reassignment to correct group",
       title = "(b)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0),
                   labels = c("Native", "St. Croix Falls Strain")) +
  theme_classic() +
  theme(plot.margin = margin(t = 5,  # Top margin
                             r = 10,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))

hybrid_plots <- density_plot + hybrid_reassignment_plot

ggsave(filename = "density_hybrids_2panel.pdf",
       plot = hybrid_plots,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hatch_int",
       height = 4,
       width = 10,
       units = "in")

ggsave(filename = "density_hybrids_2panel.png",
       plot = hybrid_plots,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hatch_int",
       height = 4,
       width = 10,
       units = "in")

#### Plot them on a map ####

# Read in necessary shape files
HUC8_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp")

HUC2_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Major_Basins/Major_Basins.shp")

WMU_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Water_Management_Units/Water_Management_Units.shp")

# Prep admixture df and lat long df for mapmixture function
HI_mapmixture <- round(prediction$posterior, 2) %>%
  as_tibble(rownames = "SampleID") %>% 
  mutate(WaterbodyName = Data_2205@pop) %>% 
  rename(Hatchery = "St.Croix_HD", Native = "Hybrid_native") %>% 
  select(WaterbodyName, SampleID, Native, Hatchery)

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
HI_map <- mapmixture(admixture_df = HI_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93.5, 
                                  xmax = -86.5, 
                                  ymin = 42, 
                                  ymax = 47.5),
                     cluster_names = c("St. Croix Falls Strain", "Native"),
                     cluster_cols = c("#D81B60", "#1E88E5"),
                     pie_size = 0.25,
                     pie_border = 0.05,
                     land_colour = "grey80",
                     sea_colour = NA,
                     arrow_size = 3,
                     arrow_position = "tr",
                     scalebar_size = 1.5,
                     scalebar_position = "tr",
                     #axis_title_size = 8, These get overwritten with theme_classic()
                     #axis_text_size = 6,
                     #plot_title = "Hatchery ancestry",
                     #plot_title_size = 10
                     ) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  labs(fill = "Reference\ngroup") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

ggsave(filename = "HI_map.pdf",
       plot = HI_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hatch_int",
       height = 5,
       width = 6,
       units = "in")

ggsave(filename = "HI_map.png",
       plot = HI_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hatch_int",
       height = 5,
       width = 6,
       units = "in")

##################################
# Plot for presentation
Melanc_willow_df <- Assignment_probs %>% 
  filter(WaterbodyName == "Willow Creek" |
         WaterbodyName == "Melancthon Creek")

Melanc_willow <- Melanc_willow_df %>% 
  ggplot(aes(SampleID, Probability, fill = factor(Cluster))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~HUC_8 + WaterbodyName, 
             switch = "x", 
             scales = "free", 
             space = "free") +
  labs(x = "Waterbody", y = "Admixture proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(values = c("#D81B60", "#1E88E5")) +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90),
        panel.grid = element_blank()) + 
  guides(fill = guide_legend(title = "Ancestry")) 

Melanc_willow <- Melanc_willow_df %>% 
  mutate(Cluster = fct_relevel(Cluster, c("Native", "Hatchery"))) %>%
  ggplot(aes(x = SampleID, y = Probability, fill = Cluster)) +
  geom_col(show.legend = TRUE) + 
  facet_nested(cols = vars(#HUC_8, 
                           WaterbodyName),
               switch = "x",
               nest_line = element_line(linewidth = 1, lineend = "round"),
               solo_line = TRUE,
               resect = unit(0.05, "in"),
               scales = "free", 
               space = "free") +
  labs(x = "", 
       y = "Assignment\nprobability",
       #title = "Great Lakes Region (HUC 2)",
       fill = "Simulated reference group") +
  scale_y_continuous(expand = c(0, 0),
                     position = "left",
                     breaks = seq(0, 1, by = 0.5)) +
  scale_fill_manual(values = c("#1E88E5", "#D81B60")) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        #strip.text.x = element_blank(),
        #strip.text.x = element_text(angle = -90,
        #                            hjust = 0),
        legend.position = "right",
        legend.direction = "vertical") +
  guides(fill = guide_legend(nrow = 1))

ggsave(filename = "Melanc_willow.png",
       plot = Melanc_willow,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 2,
       width = 6,
       units = "in")

