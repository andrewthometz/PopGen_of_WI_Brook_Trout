library(tidyverse)
library(patchwork)
library(sf)
library(RColorBrewer)
library(ggh4x)
#library(ggrepel)
#library(ggspatial)

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
#### Plot most supported number of clusters (K = 6) ####
# K = 6 STRUCTURE run
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

brewer.pal(n = 6, name = "Set1")
K6_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3")
#K6_colors <- c("#984EA3", "#FF7F00", "#4DAF4A", "#E41A1C", "yellow3", "#377EB8")

# Plot
K6_1 <- K6_longer %>% 
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
  scale_fill_manual(values = K6_colors) +
  #scale_fill_brewer(palette = "Set1") +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0))

K6_2 <- K6_longer %>% 
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
  scale_fill_manual(values = K6_colors) +
  #scale_fill_brewer(palette = "Set1") +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1))

K6_plots <- K6_1 / K6_2

ggsave(filename = "K6_str_plot.pdf",
       plot = K6_plots,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 10,
       width = 14,
       units = "in")

ggsave(filename = "K6_str_plot.png",
       plot = K6_plots,
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

# Figure out cluster for each pop
# Prep admixture df and lat long df for mapmixture function
K6_mapmixture <- K6 %>% 
  left_join(Samples_2205) %>% 
  select(WaterbodyName, SampleID, K1, K2, K3, K4, K5, K6) %>% 
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
K6_map <- mapmixture(admixture_df = K6_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93.5, 
                                  xmax = -86.5, 
                                  ymin = 42, 
                                  ymax = 47.5),
                     cluster_names = c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6"),
                     cluster_cols = K6_colors,
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
                     plot_title = "STRUCTURE (K = 6)",
                     plot_title_size = 10) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  labs(fill = "Cluster") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 6))

ggsave(filename = "K6_map.pdf",
       plot = K6_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "K6_map.png",
       plot = K6_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 5,
       width = 5,
       units = "in")
