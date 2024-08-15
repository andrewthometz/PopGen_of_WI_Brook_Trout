library(tidyverse)
library(patchwork)
library(sf)
library(RColorBrewer)

#devtools::install_github("Tom-Jenkins/mapmixture")
library(mapmixture)

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

############################
#### Plot them on a map ####
# Read in necessary shape files
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
brewer.pal(n = 3, name = "Set1")
K3_colors <- c("#E41A1C", "#377EB8", "#4DAF4A")

# Plot using mapmixture function
K3_map <- mapmixture(admixture_df = K3_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93, 
                                  xmax = -86.5, 
                                  ymin = 42.4, 
                                  ymax = 47.5),
                     cluster_names = c("1 (St. Croix Falls Strain)", "2", "3"),
                     cluster_cols = K3_colors,
                     pie_size = 0.25,
                     pie_border = 0.05,
                     #pie_opacity = 0.75,
                     land_colour = "grey80",
                     sea_colour = NA,
                     arrow = FALSE,
                     scalebar = FALSE,
                     axis_title_size = 8,
                     axis_text_size = 6,
                     plot_title = "(A)  K = 3",
                     plot_title_size = 10) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6))
        #axis.title.x = element_blank())
          #guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1),
          #                           nrow = 1)) +
          #labs(fill = "Cluster") +
          #theme(legend.title = element_text(size = 10),
          #      legend.text = element_text(size = 8),
          #      legend.position = "bottom",
          #      axis.title.x = element_blank())

# K = 6
K6 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Final_run_2205/AssProbs_CleanedUp/K6_AssProbs_CleanedUp.txt") %>% 
  select(-c(n, percent_miss)) %>% 
  mutate(C6 = as.numeric(C6)) %>% 
  rename(K1 = C4,
         K2 = C6,
         K3 = C3,
         K4 = C1,
         K5 = C5,
         K6 = C2)

K6_mapmixture <- K6 %>% 
  left_join(Samples_2205) %>% 
  select(WaterbodyName, SampleID, K1, K2, K3, K4, K5, K6) %>% 
  filter(WaterbodyName != "St. Croix Falls Strain")

# Create custom color palette
brewer.pal(n = 6, name = "Set1")
K6_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3")

# Plot using mapmixture function
K6_map <- mapmixture(admixture_df = K6_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93, 
                                  xmax = -86.5, 
                                  ymin = 42.4, 
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
                     plot_title = "(B)  K = 6",
                     plot_title_size = 10) +
  theme_classic() +
  theme(legend.position = "none",
        #axis.text = element_blank(),
        axis.text = element_text(size = 6))
        #axis.title = element_blank())
  #guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1),
  #                           nrow = 2)) +
  #labs(fill = "Cluster") +
  #theme(legend.title = element_text(size = 10),
  #      legend.text = element_text(size = 8),
  #      legend.position = "bottom",
  #      legend.direction = "horizontal",
  #      axis.text.y = element_blank(),
  #      axis.title.y = element_blank())

# K = 9
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

K9_mapmixture <- K9 %>% 
  left_join(Samples_2205) %>% 
  select(WaterbodyName, SampleID, K1, K2, K3, K4, K5, K6, K7, K8, K9) %>% 
  filter(WaterbodyName != "St. Croix Falls Strain")

# Create custom color palette
brewer.pal(n = 9, name = "Set1")
K9_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3", "#A65628", "#F781BF", "#999999")

# Plot using mapmixture function
K9_map <- mapmixture(admixture_df = K9_mapmixture, 
                     coords_df = Lats_Longs,
                     basemap = HUC2_shp,
                     boundary = c(xmin = -93, 
                                  xmax = -86.5, 
                                  ymin = 42.4, 
                                  ymax = 47.5),
                     cluster_names = c("1 (St. Croix Falls Strain)", "2", "3", "4", "5", "6", "7", "8", "9"),
                     cluster_cols = K9_colors,
                     pie_size = 0.25,
                     pie_border = 0.05,
                     #pie_opacity = 0.75,
                     land_colour = "grey80",
                     sea_colour = NA,
                     arrow = FALSE,
                     scalebar = FALSE,
                     axis_title_size = 8,
                     axis_text_size = 6,
                     plot_title = "(C)  K = 9",
                     plot_title_size = 10) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1),
                             ncol = 1)) +
  labs(fill = "Cluster\n(STRUCTURE)") +
  theme_classic() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        #legend.position = "right",
        #legend.direction = "vertical",
        axis.text = element_text(size = 6))
        #axis.text.y = element_blank(),
        #axis.title = element_blank())

# Join the 3 maps together
STRUCTURE_maps <- K3_map + K6_map + K9_map + guide_area() + plot_layout(guides = "collect")

ggsave(filename = "STRUCTURE_maps_3panel.pdf",
       plot = STRUCTURE_maps,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 8,
       width = 7,
       units = "in")

ggsave(filename = "STRUCTURE_maps_3panel.png",
       plot = STRUCTURE_maps,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Genetic_structure/STRUCTURE",
       height = 8,
       width = 7,
       units = "in")

