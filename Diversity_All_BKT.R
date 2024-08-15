library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(radiator)
library(sf)
library(ggspatial)
library(ggOceanMaps)

#### Read in Master brook trout genepop file ####
UNIFIED_BKT <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/UNIFIED_BKT_genepop.gen",
                            ncode = 3L,
                            quiet = FALSE)

nAll(UNIFIED_BKT) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Marker") %>% 
  as_tibble() %>% 
  rename("n_alleles" = ".") %>% 
  count(n_alleles) %>% 
  summarize(mean_alleles = mean(n_alleles),
            max_alleles = max(n_alleles),
            min_alleles = min(n_alleles))

#### Read in metadata ####
# 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls domestic",
         HUC_2 = "Hatchery",
         HUC_4 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_4, HUC_2)

# 2205 metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  #filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  #arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) %>% 
  mutate(Data_source = "Thometz")

# Erdman
Erdman_samples <- read_excel("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx") %>%
  mutate(Data_source = "Erdman")

Erdman_samples %>% 
  filter(WaterbodyName %in% Samples_2205$WaterbodyName) %>% 
  count(WaterbodyName)

# Bind the metadata
All_metadata <- Samples_2205 %>% 
  bind_rows(Erdman_samples) %>% 
  select(SampleID, WBIC, WaterbodyName, HUC_4, HUC_8, Data_source, Latitude, Longitude) %>% 
  filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) 

# Assign pop slot
UNIFIED_BKT@pop <- as_factor(All_metadata$WaterbodyName)

#### Calculate genetic diversity measures ####
library(hierfstat)

gd <- basic.stats(UNIFIED_BKT)

H_expected <- as.data.frame.matrix(gd$Hs) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:151,
               names_to = "pop",
               values_to = "He") %>% 
  drop_na(He) %>% 
  group_by(pop) %>% 
  summarize("He" = round(mean(He), 3))

H_observed <- as.data.frame.matrix(gd$Ho) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:151,
               names_to = "group",
               values_to = "Ho") %>% 
  drop_na(Ho) %>% 
  group_by(group) %>% 
  summarize("Ho" = round(mean(Ho), 3))

F_is <- as.data.frame.matrix(gd$Fis) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:151,
               names_to = "group",
               values_to = "Fis") %>% 
  drop_na(Fis) %>% 
  group_by(group) %>% 
  summarize("Fis" = round(mean(Fis), 3))

Allelic_richness <- as.data.frame(allelic.richness(UNIFIED_BKT)) %>% 
  as_tibble(rownames = "locus") %>% 
  select(-min.all) %>% 
  pivot_longer(2:151,
               names_to = "group",
               values_to = "ar") %>% 
  drop_na(ar) %>% 
  group_by(group) %>% 
  summarize("ar" = round(mean(ar), 2))

GD_tibble <- bind_cols(H_expected$pop,
                       Allelic_richness$ar, 
                       H_observed$Ho, 
                       H_expected$He,
                       F_is$Fis) %>% 
  rename("WaterbodyName" = 1,
         "Ar" = 2,
         "Ho" = 3, 
         "He" = 4,
         "Fis" = 5)

#write.csv(GD_tibble, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Analyses/Genetic_diversity/GD_Unified_BKT.csv")

#### Map it ####
# Read in necessary shape files
# Read in necessary shape files
HUC8_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp")

HUC2_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Major_Basins/Major_Basins.shp")

HUC4_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/NHD_H_Wisconsin_State_Shape/Shape/WBDHU4.shp") %>% 
  filter(name %in% Samples_2205$HUC_4)

HUC4_clipped <- clip_shapefile(HUC4_shp,
                               limits = HUC2_shp,
                               return.boundary = FALSE) # Select TRUE to retain metadata, must subset later with $


# Grab lats longs and convert to stat sf or whatever
Map_df <- All_metadata %>% 
  select(WaterbodyName, Latitude, Longitude) %>% 
  drop_na(Latitude, Longitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>% 
  distinct() %>% 
  right_join(GD_tibble) #%>% 
  #filter(geometry != "POINT EMPTY")

# Plot map
Ar_map <- HUC4_clipped %>% 
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
  geom_sf(data = Map_df,
          aes(geometry = geometry, color = Ar),
          size = 1.5) +
  scale_color_gradient(low = "red", 
                       high = "blue",
                       breaks = c(min(Map_df$Ar),
                                  round(mean(Map_df$Ar), 2),
                                  max(Map_df$Ar)),
                       name = bquote(A[r])) +
  annotation_scale(location = "tr") +
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.0, "in"), 
                         pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_classic()

ggsave(filename = "Ar_map_unified_BKT.pdf",
       plot = Ar_map,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 5,
       width = 5,
       units = "in")

ggsave(filename = "Ar_map_unified_BKT.png",
       plot = Ar_map,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 5,
       width = 5,
       units = "in")

#### Model lat vs Ar ####

model_df <- GD_tibble %>% 
  left_join(All_metadata) %>% 
  drop_na(Latitude)

Ar_mod_lat <- lm(Ar ~ Latitude, data = model_df)
par(mfrow = c(2,2))
plot(Ar_mod_lat)
summary(Ar_mod_lat)

Ar_lat_plot <- model_df %>% 
  ggplot(aes(x = Latitude, y = Ar)) + 
  geom_smooth(method = "glm") +
  geom_point() +
  xlab("Latitude") +
  ylab("Allelic richness") +
  scale_x_continuous(expand = c(0,0.02)) +
  scale_y_continuous(expand = c(0,0.002)) +
  annotate(geom = "text",
           label = bquote(R^2~"= 0.059"),
           x = 45.5, 
           y = 3) +
  annotate(geom = "text",
           label = bquote(p-value~"= 2.2e-16"),
           x = 45.5, 
           y = 2.8) +
  theme_classic()

ggsave(filename = "Ar_lat_regression_unified_BKT.pdf",
       plot = Ar_lat_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 5,
       width = 6,
       units = "in")

ggsave(filename = "Ar_lat_regression_unified_BKT.png",
       plot = Ar_lat_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 5,
       width = 6,
       units = "in")
