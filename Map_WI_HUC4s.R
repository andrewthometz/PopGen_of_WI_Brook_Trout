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
library(ggOceanMaps)

library(ggmap)

#### Grab lats longs from Samples_2205 ####
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  select(WaterbodyName, HUC_4, HUC_8, HUC_12, Latitude, Longitude, WBIC) %>% 
  distinct()

# Read in necessary shape files
HUC8_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp")

HUC2_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Major_Basins/Major_Basins.shp")

#WMU_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Water_Management_Units/Water_Management_Units.shp")

DriftlessArea <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/FHP_DARE_Boundary_2013.shp/FHP_DARE_Boundary_2013.shp")

DA_clipped <- clip_shapefile(DriftlessArea,
                             limits = HUC2_shp,
                             return.boundary = FALSE)

HUC4_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/NHD_H_Wisconsin_State_Shape/Shape/WBDHU4.shp") %>% 
  filter(name %in% Samples_2205$HUC_4)

HUC4_clipped <- clip_shapefile(HUC4_shp,
                               limits = HUC2_shp,
                               return.boundary = TRUE) # Select TRUE to retain metadata, must subset later with $

##############
temp_metadata <- Samples_2205 %>% 
  distinct()

my_colors <- c("#377EB8", "#4DAF4A", "#E41A1C", "#984EA3", "#A65628", "yellow3", "#F781BF", "#FF7F00")

my_breaks <- temp_metadata %>% 
  select(HUC_4) %>% 
  distinct() %>% 
  drop_na()

# Standard map
HUC_DA_map <- HUC4_clipped$shapefile %>% 
  ggplot() +
  geom_sf(aes(fill = HUC4_clipped$shapefile$name),
          #alpha = 0.5,
          color = "grey40",
          linewidth = 0.5
          ) +
  geom_sf(data = HUC2_shp,
          fill = NA,
          #alpha = 0.75,
          color = "black",
          linewidth = 1) +
  geom_sf(data = DA_clipped,
          linetype = "dashed",
          fill = NA,
          #alpha = 0.15,
          color = "black",
          linewidth = 1) +
  #geom_sf(data = GD_map_2205,
          #aes(geometry = geometry),
          #color = "black",
          #fill = "white",
          #size = 1.5) +
  #annotate(geom = "text",
  #         label = "Driftless\nArea",
  #         x = -92.1, 
  #         y = 43.6,
  #         size = 6) +
  scale_fill_manual(#name = "Subregion (HUC 4)",
                    #na.value = "black",
                    values = my_colors,
                    breaks = my_breaks$HUC_4) +
  #labs(x = "Longitude",
  #     y = "Latitude") +
  #annotation_scale(location = "tr") +
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.5, "in"), 
                         pad_y = unit(1, "in"),
                         style = north_arrow_orienteering) +
  #theme_classic() +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
  

ggsave(filename = "HUCs_Driftless_map.pdf",
       plot = HUC_DA_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "HUCs_Driftless_map.png",
       plot = HUC_DA_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 5,
       width = 5,
       units = "in")
