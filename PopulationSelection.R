library(readxl)
library(tidyverse)
library(sf)
library(lubridate)
library(timetk)
library(chron)
library(formattable)


###### MCGL Database ######

### MCGL Samples ###
colTypes <-
  c(rep("guess", 28),
    "text",
    "numeric",
    "numeric",
    "guess",
    "guess",
    "text")

## Read in individual databases
mcgl1 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Archive (04-00001 to 10-18506).xlsx",
                    col_types = colTypes,
                    range = cell_cols("A:AH")) %>% 
  select(SampleID, CommonName, Species, YearCollected, WaterbodyName, State, County, Latitude, Longitude, Comments)

mcgl2 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (11-00001 to 17-15647).xlsx",
                    col_types = colTypes,
                    range = cell_cols("A:AH")) %>% 
  select(SampleID, CommonName, Species, YearCollected, WaterbodyName, State, County, Latitude, Longitude, Comments)

mcgl3 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (18-00001 to 19-28271).xlsx",
                    col_types = colTypes,
                    range = cell_cols("A:AH")) %>% 
  select(SampleID, CommonName, Species, YearCollected, WaterbodyName, State, County, Latitude, Longitude, Comments)

mcgl4 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (21-00001 to 22-19350).xlsx",
                    col_types = colTypes,
                    range = cell_cols("A:AH")) %>% 
  select(SampleID, CommonName, Species, YearCollected, WaterbodyName, State, County, Latitude, Longitude, Comments)

# bind the four databases together
mcglBKT <- bind_rows(mcgl1, mcgl2, mcgl3, mcgl4) 

# Filter to just Wisconsin brook trout with a known stream name
mcglBKT.clean <- mcglBKT %>% 
  filter(CommonName == "Brook trout" | CommonName == "brook trout" | CommonName == "Brook Trout") %>% 
  filter(State == "WI" | State == "Wisconsin") %>% 
  add_count(WaterbodyName) %>% 
  distinct(WaterbodyName, .keep_all = TRUE) %>% 
  select(SampleID, WaterbodyName, YearCollected, County, Latitude, Longitude, Comments, n)

############################################################################################################################################
# Read in Erdman records
erdmanPops <- read_excel("C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Erdman_BKT_Genotypes.xlsx") %>% 
  select(PopCode, WaterbodyName, Data_Source, Latitude, Longitude) %>% 
  filter(Data_Source == "Erdman et al. (2020)",
         WaterbodyName != "Manchester Fish Hatchery") %>% 
  distinct(WaterbodyName, .keep_all = TRUE) %>% 
  arrange(WaterbodyName) %>% 
  mutate(Longitude = as.numeric(Longitude),
         Latitude = as.numeric(Latitude), .keep = "unused") %>% 
  st_as_sf(.,
           coords = c("Longitude", "Latitude"),
           crs = 4326)

### Print out of Erdman Pops ###

write.csv(erdmanPops, "C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Erdman_Pops.csv")

### Find BKT samples in MCGL that were not sampled in Erdman et al. 2022 ###

usable.MCGL.samples <- mcglBKT.clean %>% 
  left_join(erdmanPops, by = "WaterbodyName") %>% 
  select(SampleID, WaterbodyName, County, n, PopCode, YearCollected) %>% 
  filter(is.na(PopCode)) %>%    # this should remove most of Erdman pops as they all have abbreviations
  drop_na(County) %>% 
  filter(YearCollected != "2022") %>% 
  filter(n >= 30) %>% 
  filter(WaterbodyName != "18 Mile Creek",
         WaterbodyName != "Esofea Branch North Fork Bad Axe River",
         WaterbodyName != "Big Roche-a-Cri Creek",
         WaterbodyName != "East Branch Eau Claire River",  # manually removed these due to slight name differences
         WaterbodyName != "KC Creek",
         WaterbodyName != "Markgraf Springs",
         #WaterbodyName != "Gilbert Creek",
         WaterbodyName != "Pompey Pillar",
         WaterbodyName != "Prairie River",
         WaterbodyName != "Seven Mile Creek",
         WaterbodyName != "South Fork La Crosse River",
         WaterbodyName != "Squaw Creek",
         WaterbodyName != "Tainter Creek",
         #WaterbodyName != "St. Croix",
         #WaterbodyName != "St. Croix Falls Fish Hatchery",
         #WaterbodyName != "St. Croix River",
         WaterbodyName != "Osceola Fish Hatchery",
         #WaterbodyName != "Little Willow Creek",
         WaterbodyName != "Oconto River",
         WaterbodyName != "Tagatz Creek",                #Tagatz Creek already in list
         WaterbodyName != "Strutt Creek",
         WaterbodyName != "Woods Creek",            #Woods Creek already in list
         WaterbodyName != "Cutler Creek",
         WaterbodyName != "Hay River",
         WaterbodyName != "Clam River Flowage",
         WaterbodyName != "Kinnickinnic River",
         WaterbodyName != "Little Willow Creek",    #Little Willow already in list
         WaterbodyName != "Gilbert Creek",
         WaterbodyName != "Pine Creek",
         WaterbodyName != "Duncan Creek",
         WaterbodyName != "Shioc River",
         WaterbodyName != "Lowery Creek",
         WaterbodyName != "West Branch Mill Creek",
         WaterbodyName != "Hynek Creek",
         WaterbodyName != "Upper Pine Unnamed Tributary",
         WaterbodyName != "Story Creek") %>%
  select(-PopCode) %>% 
  arrange(County)

write.csv(usable.MCGL.samples, "C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/MCGLsamples_not_in_Erdman.csv")

##########################################################################################################################################
###### Population Selection ######

# Read in individual databases

Statewide_BKT_CPUE <- read_excel("C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Statewide brook trout CPUE.xlsx", sheet = "Sheet1") %>% 
  select("WBIC", "Survey_Year", "n_caught", "Hours", "CPE_Hour", "Station_Name") %>% 
  group_by(WBIC) %>% 
  filter(Survey_Year == max(Survey_Year)) %>% 
  filter(n_caught == max(n_caught))

Statewide_BKT_Stocking <- read_excel("C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Statewide brook trout CPUE.xlsx", sheet = "Stocking data") %>% 
  select("WBIC", "Stocking_Year", "Strain", "n_stocked", "Stock_Source_Group") %>% 
  group_by(WBIC) %>% 
  filter(Stocking_Year == max(Stocking_Year)) #%>% 
  #filter(n_stocked == max(n_stocked)) #%>% 
  #rename(Number_Caught = Number_of_Fish)

#Stocking_Hist <- read_excel("C:/Users/27tho/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/northern_counties_stocking_data_all_years.xlsx") %>% 
#  select("WBIC", "Number_Fish_Stocked", "Year", "Strain_Stock") %>% 
#  group_by(WBIC) %>% 
#  filter(Year == max(Year)) %>% 
#  rename(Most_Recent_Stocking = Year) #%>% 
  #filter(Number_Fish_Stocked == max(Number_Fish_Stocked))
                    
BKT_pop_list <- read_excel("C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Thometz.PotentialPopulationsDatabase.xlsx") %>% 
  select("WaterbodyName", "Collectors", "n_caught", "Shock_time", "Status", "Watershed", "County", "HUC8_Subbasin", "WBIC", 
         "Latitude", "Longitude", "In_DNR_Schedule", "DNR_Contact", "Warden_Contact", "Warden_Number", "Erdman", "MCGL_Samples") %>% 
  mutate(Shock_time_2 = ms(c(Shock_time)), .keep = "unused") %>% 
  mutate(Hours = as.numeric(Shock_time_2, "hours"), .keep = "unused") %>% 
  mutate(CPE_Hour = (n_caught/Hours))

# Join the databases

Joined_Database <- BKT_pop_list %>% 
  left_join(Statewide_BKT_CPUE, by = "WBIC") %>%  
  left_join(Statewide_BKT_Stocking, by = "WBIC") %>% 
  mutate(n_caught = coalesce(n_caught.x, n_caught.y),
         CPE_Hour = coalesce(CPE_Hour.x, CPE_Hour.y),
         Hours = coalesce(Hours.x, Hours.y), .keep = "unused") %>% 
  mutate(CPE_Hour = round(CPE_Hour, digits = 2),
         Hours = round(Hours, digits =2), .keep = "unused")

AvailablePops <- Joined_Database %>% 
  filter(Collectors == "Available") %>% 
  drop_na(Survey_Year)
#         CPE_Hour >= 40)

######################################################################################################################################

Final_List <- Joined_Database %>% 
  filter(Collectors == "WDNR (2022)" | Collectors == "WDNR (older MCGL samples)" | Collectors == "A. Thometz") %>% 
  distinct() %>% 
  select("WaterbodyName", "Collectors", "Watershed", "County", "WBIC", "Latitude", "Longitude", "Stocking_Year", 
         "Stock_Source_Group", "Strain", "n_stocked", "Survey_Year", "n_caught", "Hours", "CPE_Hour") %>%
  #mutate(CPE_Hour = case_when(is.na(CPE_Hour) ~ "No data", TRUE ~ as.character(CPE_Hour)))
  st_as_sf(.,
           coords = c("Longitude", "Latitude"),
           crs = 4326)

write.csv(Final_List, "C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Updated_BKT_Survey.csv")

Lab_list <- Final_List %>% 
  left_join(mcglBKT.clean, by = "WaterbodyName") %>% 
  select("WaterbodyName", "WBIC", "SampleID") %>% 
  arrange(SampleID)

write.csv(Lab_list, "C:/Users/athom393/OneDrive - UWSP/Project data/bktSurvey Ch. 2/PopulationSelection/Lab_list.csv")


#################################################################################################################################### 
###### Show candidates on map ######

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(ggspatial)
library(cowplot)
library(paletteer)

# Download file, extract all, copy text address of file, tab down to .shp file
rivers_streams_SHP <- read_sf("C:/Users/athom393/OneDrive - UWSP/R_Thesis/24k_Hydro_Flowlines_(Rivers_Streams)/24k_Hydro_Flowlines_(Rivers_Streams).shp") %>% 
  st_transform(crs = 4326)

Rivers_Streams <- rivers_streams_SHP %>% 
  filter(ROW_NAME != "Unnamed",
         STREAM_ORD >= 4)

HUC8_SHP <- read_sf("C:/Users/athom393/OneDrive - UWSP/R_Thesis/Hydrologic_Units_-_8_digit_(Subbasins)/Hydrologic_Units_-_8_digit_(Subbasins).shp") %>% 
  st_transform(crs = 4326)

HUC8 <- HUC8_SHP %>% 
  filter(HUC8_NAME != "Lake Superior",
         HUC8_NAME != "Lake Michigan",
         HUC8_NAME != "St. Louis",
         HUC8_NAME != "Lower Rock",
         HUC8_NAME != "Des Plaines",
         HUC8_NAME != "Pike-Root",
         HUC8_NAME != "Apple-Plum",
         HUC8_NAME != "Kishwaukee")

# Get Wisconsin state data
install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")


states <- ne_states(country = "United States of America",
                    returnclass = "sf")

WI <- states %>%   
  filter(name == "Wisconsin")

#combine my survey list with Erdman pops for map

Mine_plus_Erdman <- Final_List %>% bind_rows(erdmanPops)

# Run a pipeline to make the map
WI %>% 
  ggplot() +
  geom_sf(data = WI,
          fill = "azure3") +
  geom_sf(data = Rivers_Streams,
          color = "azure4",
          alpha = 0.4,
          size = 0.1) +
  geom_sf(data = HUC8, 
          alpha = 0.01, 
          color = "white", 
          size = 0.5) +
  geom_sf(data = Final_List,
          aes(geometry = geometry),
          color = "blue",
          size = 1) +
  geom_sf(data = erdmanPops,
          aes(geometry = geometry),
          color = "red",
          size = 1) +
  #geom_sf_text(data = Final_List,
               #mapping = aes(label = Final_List$Population),
               #size = 2) +
  #scale_color_manual(labels = c("MCGL samples", "Sample with Mitro", "Sample with own crew", "Sample by DNR"), 
  #                   values = c("Blue", "Yellow", "Red", "springgreen3"),
  #                   name = "Sample Acquisition") +
  annotation_scale(location = "tr") +
  annotation_north_arrow(location = "bl",
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_orienteering()) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Brook trout survey populations") +
  labs(caption = "") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.35), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot")
