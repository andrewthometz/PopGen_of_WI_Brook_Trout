library(tidyverse)
library(patchwork)
library(sf)
#library(ggrepel)
#library(ggspatial)
library(RColorBrewer)
library(ggh4x)

### For Evanno method ###
# Clear environment and restart session before running Pophelper
#library(pophelper)
#library(pophelperShiny)
#runPophelper()

#devtools::install_github("Tom-Jenkins/mapmixture")
library(mapmixture)
#launch_mapmixture()

# Prep 2111 data to work with plotting
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain",
         HUC_2 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_2)

# Read in 2205 metadata data
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_2, Latitude, Longitude) %>% 
  bind_rows(Samples_2111)

#########################################################
#### Plot most supported number of clusters (K = 9) ####
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

brewer.pal(n = 9, name = "Set1")
K9_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3", "#A65628", "#F781BF", "#999999")

# Plot
K9_1 <- K9_longer %>% 
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
  scale_fill_manual(values = K9_colors) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0))

K9_2 <- K9_longer %>% 
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
       fill = "Cluster (STRUCTURE)") +
  scale_y_continuous(expand = c(0, 0),
                     position = "left",
                     breaks = seq(0, 1, by = 0.5)) +
  scale_fill_manual(values = K9_colors) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1))

K9_plots <- K9_1 / K9_2

ggsave(filename = "K9_str_plot.pdf",
       plot = K9_plots,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 10,
       width = 14,
       units = "in")

ggsave(filename = "K9_str_plot.png",
       plot = K9_plots,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
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
K9_mapmixture <- K9 %>% 
  left_join(Samples_2205) %>% 
  select(WaterbodyName, SampleID, K1, K2, K3, K4, K5, K6, K7, K8, K9) %>% 
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
K9_map <- mapmixture(admixture_df = K9_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93.5, 
                                  xmax = -86.5, 
                                  ymin = 42, 
                                  ymax = 47.5),
                     cluster_names = c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7", "8", "9"),
                     cluster_cols = K9_colors,
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
                     plot_title = "STRUCTURE (K = 9)",
                     plot_title_size = 10) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(fill = "Cluster") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 6))

ggsave(filename = "K9_map.pdf",
       plot = K9_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "K9_map.png",
       plot = K9_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 5,
       width = 5,
       units = "in")
