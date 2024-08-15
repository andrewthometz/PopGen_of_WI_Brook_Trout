library(tidyverse)
library(patchwork)
library(sf)
library(RColorBrewer)
#library(ggrepel)
#library(ggspatial)
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
#### Plot most supported number of clusters (K = 3) ####
# K = 3 STRUCTURE run
K3 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Final_run_2205/AssProbs_CleanedUp/K3_AssProbs_CleanedUp.txt") %>% 
  select(-1) %>% 
  mutate(C3 = as.numeric(C3)) %>% 
  pivot_longer(cols = C1:C3, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205) %>% 
  mutate(Cluster = str_replace(Cluster, "C", ""),
         Cluster = str_replace(Cluster, "1", "1 (St. Croix Falls Strain)"),
         Cluster = fct_relevel(Cluster, c("1 (St. Croix Falls Strain)", "2", "3")),
         HUC_8 = fct_relevel(HUC_8, "Hatchery", after = Inf))

# Plot  
K3_1 <- K3 %>% 
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
  #scale_fill_manual(values = rainbow(3)) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0))

K3_2 <- K3 %>% 
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
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() + 
  theme(panel.spacing = unit(0.1, "line"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = -90,
                                    hjust = 0),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1))

K3_plots <- K3_1 / K3_2

ggsave(filename = "K3_str_plot.pdf",
       plot = K3_plots,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 10,
       width = 14,
       units = "in")

ggsave(filename = "K3_str_plot.png",
       plot = K3_plots,
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
K3_mapmixture <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Final_run_2205/AssProbs_CleanedUp/K3_AssProbs_CleanedUp.txt") %>% 
  select(-1) %>% 
  left_join(Samples_2205) %>% 
  select(WaterbodyName, SampleID, C1, C2, C3) %>% 
  mutate(C3 = as.numeric(C3)) %>% 
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

# Create custom color palette
Set1_colors <- brewer.pal(n = 3, name = "Set1")

# Plot using mapmixture function
K3_map <- mapmixture(admixture_df = K3_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93.5, 
                                  xmax = -86.5, 
                                  ymin = 42, 
                                  ymax = 47.5),
                     cluster_names = c("1 (St. Croix Falls Strain)", "2", "3"),
                     cluster_cols = Set1_colors,
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
                     plot_title = "STRUCTURE (K = 3)",
                     plot_title_size = 10) +
          guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
          labs(fill = "Cluster") +
          theme(legend.title = element_text(size = 10),
                legend.text = element_text(size = 8))

ggsave(filename = "K3_map.pdf",
       plot = K3_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "K3_map.png",
       plot = K3_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 5,
       width = 5,
       units = "in")
