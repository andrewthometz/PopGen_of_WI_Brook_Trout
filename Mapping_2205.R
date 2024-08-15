library(tidyverse)
library(readxl)
library(adegenet)

#### Mapping ####
library(sf)
library(ggrepel)
library(ggspatial)
library(cowplot)
library(paletteer)
library(patchwork)

library(ggmap)

#### Grab lats longs from Samples_2205 ####
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  select(WaterbodyName, HUC_8, HUC_4, HUC_12, Latitude, Longitude, WBIC) %>% 
  distinct()

# Read in master gd file
GD_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Genetic_diversity/Diversity_2205.csv") %>% 
  left_join(Samples_2205)

# Read in necessary shape files
HUC8_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp")

HUC2_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Major_Basins/Major_Basins.shp")

HUC4_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/NHD_H_Wisconsin_State_Shape/Shape/WBDHU4.shp") %>% 
  filter(name %in% Samples_2205$HUC_4)

HUC4_clipped <- clip_shapefile(HUC4_shp,
                               limits = HUC2_shp,
                               return.boundary = FALSE) # Select TRUE to retain metadata, must subset later with $
# Prep map data
GD_map_2205 <- GD_2205 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>% 
  mutate(Ar = round(Ar, digits = 2),
         Ho = round(Ho, digits = 2),
         He = round(He, digits = 2),
         Fis = round(Fis, digits = 2))

##############
# Standard map
Standard_map <- WMU_shp %>% 
  ggplot() +
  geom_sf(fill = NA,
          alpha = 0.5,
          color = "grey",
          linewidth = 0.25) +
  geom_sf(data = HUC2_shp,
          fill = NA,
          alpha = 0.75,
          color = "grey50",
          linewidth = 0.5) +
  geom_sf(data = GD_map_2205,
          aes(geometry = geometry),
          #color = "black",
          #fill = "white",
          size = 1.5) +
  labs(x = "Longitude",
       y = "Latitude") +
  annotation_scale(location = "tr") +
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.0, "in"), 
                         pad_y = unit(0.4, "in"),
                         style = north_arrow_orienteering) +
  theme_classic()

ggsave(filename = "Standard_sites_map.pdf",
       plot = Standard_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/GD_maps",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "Standard_sites_map.png",
       plot = Standard_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/GD_maps",
       height = 5,
       width = 5,
       units = "in")

###################
# Four-panel GD map
# Allelic richness
Ar_map <- HUC2_shp %>% 
  ggplot() +
  geom_sf(fill = NA,
          alpha = 0.5,
          #color = "grey50",
          linewidth = 0.25) +
  geom_sf(data = GD_map_2205,
          aes(geometry = geometry, color = Ar),
          size = 1.5) +
  scale_color_gradient(low = "red", 
                       high = "blue",
                       breaks = c(min(GD_map_2205$Ar),
                                  round(mean(GD_map_2205$Ar), 2),
                                  max(GD_map_2205$Ar)),
                       name = bquote(A[r])) +
  annotation_scale(location = "tr") +
  #annotation_north_arrow(location = "tr", 
  #                       which_north = "true", 
  #                       pad_x = unit(0.0, "in"), 
  #                       pad_y = unit(0.5, "in"),
  #                       style = north_arrow_fancy_orienteering) +
  labs(#x = "Longitude",
       y = "Latitude",
       title = "(A)") +
  #guides(color = guide_legend(title = "Allelic\nrichness",
  #                            reverse = TRUE)) +
  theme_classic() +
  theme(axis.text.x = element_blank())

# Inbreeding coefficient (Fis) (0 = low inbreeding, 1 = high inbreeding)
Fis_map <- HUC2_shp %>% 
  ggplot() +
  geom_sf(fill = NA,
          alpha = 0.5,
          #color = "grey50",
          linewidth = 0.25) +
  geom_sf(data = GD_map_2205,
          aes(geometry = geometry, color = Fis),
          size = 1.5) +
  scale_color_gradient(low = "blue", 
                       high = "red",
                       breaks = c(min(GD_map_2205$Fis),
                                  round(mean(GD_map_2205$Fis), 2),
                                  max(GD_map_2205$Fis)),
                       name = bquote(F[IS])) +
  #annotation_scale(location = "tr") +
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.0, "in"), 
                         pad_y = unit(0.05, "in"),
                         style = north_arrow_orienteering) +
  labs(#x = "Longitude",
       #y = "Latitude",
       title = "(B)") +
  #guides(color = guide_legend(title = "Inbreeding\ncoefficient",
  #                            reverse = TRUE)) +
  theme_classic() +
  theme(axis.text = element_blank())

# Observed hetereozygosity
Ho_map <- HUC2_shp %>% 
  ggplot() +
  geom_sf(fill = NA,
          alpha = 0.5,
          #color = "grey50",
          linewidth = 0.25) +
  geom_sf(data = GD_map_2205,
          aes(geometry = geometry, color = Ho),
          size = 1.5) +
  scale_color_gradient(low = "red", 
                       high = "blue",
                       breaks = c(min(GD_map_2205$Ho),
                                  round(mean(GD_map_2205$Ho), 2),
                                  max(GD_map_2205$Ho)),
                       name = bquote(H[o])) +
  #annotation_scale(location = "tr") +
  #annotation_north_arrow(location = "tr", 
  #                       which_north = "true", 
  #                       pad_x = unit(0.0, "in"), 
  #                       pad_y = unit(0.5, "in"),
  #                       style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude",
       y = "Latitude",
       title = "(C)") +
  #guides(color = guide_legend(title = "Observed\nheterozygosity",
  #                            reverse = FALSE)) +
  theme_classic() +
  theme(legend.key.size = unit(0.75, "cm"))

# Expected hetereozygosity
He_map <- HUC2_shp %>% 
  ggplot() +
  geom_sf(fill = NA,
          alpha = 0.5,
          #color = "grey50",
          linewidth = 0.25) +
  geom_sf(data = GD_map_2205,
          aes(geometry = geometry, color = He),
          size = 1.5) +
  scale_color_gradient(low = "red", 
                       high = "blue",
                       breaks = c(min(GD_map_2205$He), 
                                  round(mean(GD_map_2205$He), 2),
                                  max(GD_map_2205$He)),
                       name = bquote(H[e])) +
  #annotation_scale(location = "tr") +
  #annotation_north_arrow(location = "tr", 
  #                       which_north = "true", 
  #                       pad_x = unit(0.0, "in"), 
  #                       pad_y = unit(0.5, "in"),
  #                       style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude",
       title = "(D)") +
       #y = "Latitude") +
  #guides(color = guide_legend(title = "Expected\nheterozygosity",
  #                            reverse = FALSE)) +
  theme_classic() +
  theme(axis.text.y = element_blank())

GD_maps <- (Ar_map + Fis_map) / (Ho_map + He_map)

ggsave(filename = "GD_maps_4panel.pdf",
       plot = GD_maps,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/GD_maps",
       height = 8,
       width = 8,
       units = "in")

ggsave(filename = "GD_maps_4panel.png",
       plot = GD_maps,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/GD_maps",
       height = 8,
       width = 8,
       units = "in")

########
# Ne map (middle value is median not mean)
samples_temp <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  select(SampleID, WBIC)

Ne_estimates <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Genetic_diversity/Ne/Ne_All_63pops_LD_CleanedUp.txt") %>% 
  select(SampleID, Ne) %>% 
  left_join(samples_temp) %>% 
  select(-SampleID) %>% 
  mutate(Ne_estimate = case_when(Ne < 0 ~ "Inf",
                                 .default = as.character(Ne)),
         .keep = "unused")
  

Ne_df <- GD_map_2205 %>% 
  left_join(Ne_estimates) %>% 
  filter(Ne_estimate != "Inf") %>%
  mutate(Ne_estimate = as.numeric(Ne_estimate),
         Ne_estimate = case_when(#Ne == "Inf" ~ 1000,
                                 Ne_estimate > 500 ~ 500,
                                 .default = Ne_estimate),
         Ne_estimate = round(Ne_estimate, 0))

Ne_map <- HUC4_clipped %>% 
  ggplot() +
  geom_sf(fill = NA,
          alpha = 0.5,
          color = "grey",
          linewidth = 0.25) +
  geom_sf(data = HUC2_shp,
          fill = NA,
          alpha = 0.75,
          color = "grey50",
          linewidth = 0.5) +
  geom_sf(data = Ne_df,
          aes(geometry = geometry, color = Ne_estimate),
          size = 1.5) +
  scale_color_gradient(low = "red", 
                       high = "blue",
                       breaks = c(min(Ne_df$Ne_estimate),
                                  median(Ne_df$Ne_estimate),
                                  max(Ne_df$Ne_estimate)),
                       name = bquote(N[e])) +
  annotation_scale(location = "tr") +
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.0, "in"), 
                         pad_y = unit(0.4, "in"),
                         style = north_arrow_orienteering) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_classic()

ggsave(filename = "Ne_map.pdf",
       plot = Ne_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/GD_maps",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "Ne_map.png",
       plot = Ne_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/GD_maps",
       height = 5,
       width = 5,
       units = "in")
