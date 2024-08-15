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
Erdman_excel <- read_excel("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx") %>% 
  select(SampleID, WBIC, Pop, WaterbodyName, Latitude, Longitude, Data_source)

#### Add HUC data ####
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

# Grab HUC12 names because this file is missing WBIC's for ponds/lakes
HUC12_names <- HUC12_WBICs %>% 
  select(HUC_12, HUC12_CODE) %>% 
  distinct()

HUC12_names %>% 
  filter(HUC12_CODE == "070500030201")

# Make final df with HUC data
Erdman_data <- Erdman_excel %>% 
  mutate(HUC_2 = getHUCName(longitude = Erdman_excel$Longitude, latitude = Erdman_excel$Latitude, dataset = huc2Dat),
         HUC_4 = getHUCName(longitude = Erdman_excel$Longitude, latitude = Erdman_excel$Latitude, dataset = huc4Dat),
         HUC_6 = getHUCName(longitude = Erdman_excel$Longitude, latitude = Erdman_excel$Latitude, dataset = huc6Dat),
         HUC_8 = getHUCName(longitude = Erdman_excel$Longitude, latitude = Erdman_excel$Latitude, dataset = huc8Dat),
         HUC_10 = getHUCName(longitude = Erdman_excel$Longitude, latitude = Erdman_excel$Latitude, dataset = huc10Dat)) %>% 
  select(SampleID, WBIC, WaterbodyName, Pop, Latitude, Longitude, HUC_2, HUC_4, HUC_6, HUC_8, HUC_10) %>% 
  left_join(HUC12_WBICs, relationship = "many-to-many") %>% 
  mutate(HUC12_CODE = case_when(WaterbodyName == "Foulds Springs" ~ "070500030201",
                                WaterbodyName == "Hogelee Spring Pond 1" ~ "040302020402",
                                WaterbodyName == "Hogelee Spring Pond 2" ~ "040302020402",
                                WaterbodyName == "Krause Springs" ~ "040302020401",
                                WaterbodyName == "McGee Lake" ~ "040302020402",
                                WaterbodyName == "Rabe Lake" ~ "040302020401",
                                .default = HUC12_CODE)) %>% 
  mutate(correct_HUC12 = case_when(WaterbodyName == "Beaver Brook" ~ "070300010401",
                                   WaterbodyName == "Big Roche a Cri Creek" ~ "070700030804",
                                   WaterbodyName == "Bois Brule River" ~ "040103010705",
                                   WaterbodyName == "Chipmunk Coulee Creek" ~ "070600010502",
                                   WaterbodyName == "Comet Creek" ~ "040302021503",
                                   WaterbodyName == "East Branch Eau Claire River" ~ "070700021203",
                                   WaterbodyName == "Eighteen Mile Creek" ~ "070500070708",
                                   WaterbodyName == "Elk Creek" ~ "070400050304",
                                   WaterbodyName == "North Fork Bad Axe River" ~ "070600010305",
                                   WaterbodyName == "Gran Grae Creek" ~ "070700051803",
                                   WaterbodyName == "Little Plover River" ~ "070700030303",
                                   WaterbodyName == "Little Wolf River" ~ "040302021705",
                                   WaterbodyName == "North Branch Embarrass River" ~ "040302021202",
                                   WaterbodyName == "North Fork Clam River" ~ "070300010805",
                                   WaterbodyName == "Prairie River" ~ "070700020306",
                                   WaterbodyName == "Sawyer Creek" ~ "070300010403",
                                   WaterbodyName == "Sevenmile Creek" ~ "070700030703",
                                   WaterbodyName == "South Branch Oconto River" ~ "040301040102",
                                   WaterbodyName == "South Fork Hay River" ~ "070500070506",
                                   WaterbodyName == "South Fork La Crosse River" ~ "070400060202",
                                   WaterbodyName == "Upper Duncan Creek" ~ "070500050404",
                                   WaterbodyName == "Wausaukee River" ~ "040301080905",
                                   WaterbodyName == "Wilson Creek" ~ "070500071002",
                                  .default = HUC12_CODE)) %>%
    filter(HUC12_CODE == correct_HUC12) %>%
    left_join(HUC12_names, by = "HUC12_CODE") %>% 
    select(-c(HUC12_CODE, correct_HUC12, HUC_12.x)) %>% 
    rename(HUC_12 = HUC_12.y)

Erdman_data %>% 
  filter(is.na(HUC_12))

write_csv(Erdman_data, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_HUC_data.csv")

#### Create reusable stocking data csv ####
Stocking_data <- read_excel("X:/2205_BKT_feral_broodstock_ID/Statewide BKT info/WI_BKT_Stocking&CPUE.xlsx", 
                            sheet = "Stocking data") %>% 
  select(WBIC, WaterbodyName, Stocking_Year, n_stocked, Strain, Stock_Source_Group) %>% 
  rename(Name = WaterbodyName) %>% 
  left_join(HUC12_WBICs, relationship = "many-to-many") %>% 
  filter(HUC_12 %in% Erdman_data$HUC_12)

write_csv(Stocking_data, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_stocking_histories.csv")

#### Create reusable CPE data csv ####
CPE_data <- read_excel("X:/2205_BKT_feral_broodstock_ID/Statewide BKT info/WI_BKT_Stocking&CPUE.xlsx", 
                       sheet = ">200mile") %>% 
  select(WBIC, Survey_Year, Hours, N_Fish, CPE_Hour, CPE_Mile) %>% 
  rename(n_caught = N_Fish) %>% 
  filter(WBIC %in% Erdman_data$WBIC)

write_csv(CPE_data, "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_CPE_data.csv")
