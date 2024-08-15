library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)
library(PopGenReport)
library(radiator)

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
#HUCs <- Samples_2205 %>% 
#  select(WBIC, WaterbodyName, HUC_2, HUC_4, HUC_6, HUC_8) %>% 
#  distinct()

# Read in St.Croix hybrid data
St.Croix_HD <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/St.Croix_HD_genepop.gen", 
                            ncode = 3L, 
                            quiet = FALSE)

temp <- tibble(WaterbodyName = rep("St.Croix_HD", 100))
St.Croix_HD@pop <- as_factor(temp$WaterbodyName)

# Join them together
Joined_genind <- repool(Data_2205, St.Croix_HD)

#### Calculate individual level genetic distance ####
# Takes ~ 5 min to run
#gd.smouse(Joined_genind, verbose = TRUE) %>% 
#  tidy() %>% 
#  write.csv("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/gd.smouse.csv")

result_tidy <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/gd.smouse.csv") %>% 
  select(-1)

### Filter to find HUC12s with no bkt stocking history (10 pops) #### Needed to be manually reviewed
# Read in WDNR stocking database #
Stocking_data <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Stocking_histories_2205.csv")

never_stocked <- Samples_2205 %>% 
  filter(!(HUC_12 %in% Stocking_data$HUC_12),
         WaterbodyName != "Allequash Creek", # Allequash Springs heavily stocked
         WaterbodyName != "Bearskin Creek", # Bearskin Creek adjacent to Swamp Creek, heavily stocked
         WaterbodyName != "Becky Creek", # Becky Creek adjacent to Devils Creek, heavily stocked
         WaterbodyName != "Mt. Pelee Creek", # Mt Pelee Creek adjacent to Niebauer Springs, heavily stocked
         WaterbodyName != "Tomorrow River", # Nelsonville pond stocked several times in 70s
         WaterbodyName != "Upper Pine River", # Wild Rose Pond stocked in 1972
         WaterbodyName != "Woods Creek",
         WaterbodyName != "St. Croix Falls domestic") %>% # Woods Creek adjacent to Popple and Pine Rivers, heavily stocked
  select(SampleID, WaterbodyName, WBIC)

never_stocked %>% count(WaterbodyName)

#### Filter down to fish with greatest genetic distance from St.Croix domestic hybrids (upper 25%) ####
Distinct_bkt <- result_tidy %>% 
  filter(!str_detect(item1, "HD"),
         str_detect(item2, "HD")) %>% 
  rename(SampleID = item1,
         Domestics = item2) %>% 
  filter(SampleID %in% never_stocked$SampleID) %>% 
  left_join(never_stocked, by = "SampleID") %>% 
  group_by(SampleID) %>% 
  filter(distance == max(distance)) %>% 
  distinct(SampleID, .keep_all = TRUE) %>%
  ungroup() %>% 
  filter(distance > max(distance)*0.75) %>% # Top 25% retains reasonable quantity of fish to hybridize (209)
  mutate(Group = sample(1:2, n(), replace = TRUE)) # Randomly divide into two groups

Distinct_bkt %>% 
  group_by(WaterbodyName) %>%
  count(WaterbodyName)

Distinct_bkt %>% 
  count(SampleID)

Distinct_bkt %>% 
  count(Group)

# Bring randomized groupings back to Samples_2205
Samples_2205 <- Distinct_bkt %>% 
  select(SampleID, Group) %>% 
  right_join(Samples_2205)
  
# Subset by randomized groupings
Data_2205@pop <- as_factor(Samples_2205$Group)

Group_1 <- popsub(Data_2205, "1")

Group_2 <- popsub(Data_2205, "2")

#### Hybridize wild fish to create the population of "Native WI" brook trout ####
hybridize(Group_1, Group_2,
          n = 100,
          pop = "Hybrid_native",
          res.type = "genind",
          hyb.label = "HN") %>% 
  tidy_genind() %>% 
  write_genepop(genepop.header = "Hybridized native fish (n = 100)",
                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/Hybridized_natives")

# Read in Hybrid Native genetic data
Hybrid_natives <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Hybridized_BKT/Hybridized_natives.gen",
                               ncode = 3L, 
                               quiet = FALSE)

temp <- tibble(Hybrid_type = rep("Hybrid_native", 100))
Hybrid_natives@pop <- as_factor(temp$Hybrid_type)

#### Map the unstocked populations to ensure even coverage of state ####
library(sf)
library(ggspatial)
# Read in necessary shape files
HUC8_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp")

HUC2_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Major_Basins/Major_Basins.shp")

WMU_shp <- read_sf("X:/2205_BKT_feral_broodstock_ID/Mapping_shapefiles/Water_Management_Units/Water_Management_Units.shp")

unstocked_map <- Samples_2205 %>% 
  filter(WaterbodyName %in% never_stocked$WaterbodyName) %>% 
  select(WaterbodyName, Latitude, Longitude) %>% 
  distinct() %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326)
  
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
  geom_sf(data = unstocked_map,
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
                         style = north_arrow_fancy_orienteering) +
  theme_classic()
