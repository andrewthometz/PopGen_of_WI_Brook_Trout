library(tidyverse)
library(readxl)
library(sf)
library(adegenet)
library(vcfR)
library(ggrepel)
library(remotes)
library(sp)
#remove.packages("MazamaSpatialUtils")

install_version("MazamaSpatialUtils", version = "0.7.6", repos = "http://cran.us.r-project.org")
library(MazamaSpatialUtils)

#### Reading in MCGL Data ####
MCGL2 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (11-00001 to 17-15647).xlsx", 
                    guess_max = 1000000) %>% 
  select(SampleID, PlateID, ProjectID, CommonName, YearCollected, WaterbodyName, State, County, Collectors, Latitude, Longitude, Comments)
  
MCGL3 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (18-00001 to 19-28271).xlsx",
                    guess_max = 1000000) %>% 
  select(SampleID, PlateID, ProjectID, CommonName, YearCollected, WaterbodyName, State, County, Collectors, Latitude, Longitude, Comments)

MCGL4 <- read_excel("Z:/MCGL Database/Working/MCGL Sample Database (21-00001 to 22-19350).xlsx",
                    guess_max = 1000000) %>% 
  select(SampleID, PlateID, ProjectID, CommonName, YearCollected, WaterbodyName, State, County, Collectors, Latitude, Longitude, Comments) %>% 
  mutate(YearCollected = as.numeric(YearCollected))

MCGL_joined <- bind_rows(MCGL2, MCGL3, MCGL4) 

#### Add HUC data ####
extra_4_pops <- MCGL_joined %>% 
  filter(WaterbodyName == "Lowery Creek" |
         WaterbodyName == "Cady Creek" |
         WaterbodyName == "Melancthon Creek" |
         WaterbodyName == "South Fork Hay River") %>% 
  filter(str_detect(Comments, "WBIC")) %>% 
  mutate(Latitude = round(as.numeric(Latitude), digits = 3),
         Longitude = round(as.numeric(Longitude), digits = 3))

temporary <- MCGL_joined %>% 
  filter(str_detect(PlateID, "2205"),
         WaterbodyName != "Strutt Creek") %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>% 
  bind_rows(extra_4_pops)

setSpatialDataDir("X:/recurringProjects/spatialDat")  # Do not change
huc2Dat <- loadSpatialData('WBDHU2')
huc4Dat <- loadSpatialData('WBDHU4')
huc6Dat <- loadSpatialData('WBDHU6')
huc8Dat <- loadSpatialData('WBDHU8')
huc10Dat <- loadSpatialData('WBDHU10')

# HUC 12 data
HUC12_WBICs <- read_delim("X:/2205_BKT_feral_broodstock_ID/Statewide BKT info/HUC12_WBIC.txt") %>% 
  rename(WBIC = RIVER_SY_1,
         HUC_12 = HUC12_NAME) %>% 
  filter(WBIC != 0) %>% 
  select(WBIC, HUC_12, HUC12_CODE) %>% 
  distinct()

# Make final df with huc data
Samples_2205 <- temporary %>% 
  mutate(HUC_2 = getHUCName(longitude = temporary$Longitude, latitude = temporary$Latitude, dataset = huc2Dat),
         HUC_4 = getHUCName(longitude = temporary$Longitude, latitude = temporary$Latitude, dataset = huc4Dat),
         HUC_6 = getHUCName(longitude = temporary$Longitude, latitude = temporary$Latitude, dataset = huc6Dat),
         HUC_8 = getHUCName(longitude = temporary$Longitude, latitude = temporary$Latitude, dataset = huc8Dat),
         HUC_10 = getHUCName(longitude = temporary$Longitude, latitude = temporary$Latitude, dataset = huc10Dat)) %>% 
  mutate(WBIC = str_sub(Comments, start = 6L, end = 12L),
         WBIC = as.numeric(str_remove_all(WBIC, ","))) %>% 
  select(SampleID, WBIC, WaterbodyName, County, Latitude, Longitude, HUC_2, HUC_4, HUC_6, HUC_8, HUC_10) %>% 
  left_join(HUC12_WBICs, relationship = "many-to-many") %>% 
  mutate(correct_HUC12 = case_when(WaterbodyName == "Buena Vista Creek" ~ "070700030401",
                                   WaterbodyName == "Devils Creek" ~ "040103020304",
                                   WaterbodyName == "Elvoy Creek" ~ "040301060302",
                                   WaterbodyName == "Evergreen River" ~ "040302020303",
                                   WaterbodyName == "Flume Creek" ~ "040302021502",
                                   WaterbodyName == "Knapp Creek" ~ "070700051504",
                                   WaterbodyName == "Little Deerskin River" ~ "070700010103",
                                   WaterbodyName == "Plover River" ~ "070700030101",
                                   WaterbodyName == "Spring Brook" ~ "070700021103",
                                   WaterbodyName == "Tagatz Creek" ~ "040302010301",
                                   WaterbodyName == "Tomorrow River" ~ "040302021801",
                                   WaterbodyName == "Upper Pine River" ~ "040302022001",
                                   WaterbodyName == "Willow Creek" ~ "070700051001",
                                   WaterbodyName == "South Fork Hay River" ~ "070500070502",
                                   .default = HUC12_CODE)) %>% 
  filter(HUC12_CODE == correct_HUC12) %>% 
  select(-c(HUC12_CODE, correct_HUC12))

write_csv(Samples_2205, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv")

#### Create reusable stocking data csv ####
Stocking_data <- read_excel("X:/2205_BKT_feral_broodstock_ID/Statewide BKT info/WI_BKT_Stocking&CPUE.xlsx", 
                            sheet = "Stocking data") %>% 
  select(WBIC, WaterbodyName, Stocking_Year, n_stocked, Strain, Stock_Source_Group) %>% 
  rename(Name = WaterbodyName) %>% 
  left_join(HUC12_WBICs, relationship = "many-to-many") %>% 
  filter(HUC_12 %in% Samples_2205$HUC_12)

write_csv(Stocking_data, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Stocking_histories_2205.csv")

#### Create reusable CPE data csv ####
CPE_data <- read_excel("X:/2205_BKT_feral_broodstock_ID/Statewide BKT info/WI_BKT_Stocking&CPUE.xlsx", 
                       sheet = ">200mile") %>% 
  select(WBIC, Survey_Year, Hours, N_Fish, CPE_Hour) %>% 
  rename(n_caught = N_Fish) %>% 
  filter(WBIC %in% Samples_2205$WBIC)

write_csv(CPE_data, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CPE_data_2205.csv")
