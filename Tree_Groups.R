library(tidyverse)
library(readxl)
library(adegenet)

# Prep 2111 data to work with plotting
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain",
         HUC_2 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_2)

# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

# Add Grouping for every population
Tree_groupings <- Samples_2205 %>% 
  mutate(Tree_group = case_when(WaterbodyName == "St. Croix Falls Strain" ~ "A",
                                WaterbodyName == "Venison Creek" ~ "A",
                                WaterbodyName == "Knapp Creek" ~ "A",
                                WaterbodyName == "Little Willow Creek" ~ "A",
                                WaterbodyName == "Willow Creek" ~ "A",
                                WaterbodyName == "Allequash Creek" ~ "A",
                                WaterbodyName == "Bearskin Creek" ~ "A",
                                WaterbodyName == "Upper Pine River" ~ "A",
                                WaterbodyName == "Unnamed trib to Dell Creek (a)" ~ "A",
                                WaterbodyName == "Unnamed trib to Dell Creek (b)" ~ "A",
                                WaterbodyName == "Chambers Creek" ~ "A",
                                WaterbodyName == "Bergen Creek" ~ "A",
                                WaterbodyName == "Trout Creek" ~ "A",
                                WaterbodyName == "Chase Creek" ~ "A",
                                WaterbodyName == "Marshall Creek" ~ "A",
                                WaterbodyName == "Marshall Creek - West Branch" ~ "A",
                                WaterbodyName == "Little Deerskin River" ~ "A",
                                WaterbodyName == "Noisy Creek" ~ "A",
                                WaterbodyName == "Lepage Creek" ~ "B",
                                WaterbodyName == "Evergreen River" ~ "B",
                                WaterbodyName == "Elvoy Creek" ~ "B",
                                WaterbodyName == "Wisconsin Creek" ~ "B",
                                WaterbodyName == "Alvin Creek" ~ "B",
                                WaterbodyName == "Woods Creek" ~ "B",
                                WaterbodyName == "Swan Creek" ~ "B",
                                WaterbodyName == "Maple Creek" ~ "B",
                                WaterbodyName == "Becky Creek" ~ "B",
                                WaterbodyName == "North Otter Creek" ~ "B",
                                WaterbodyName == "Mt. Pelee Creek" ~ "B",
                                WaterbodyName == "Cap Creek" ~ "B",
                                WaterbodyName == "Spring Creek" ~ "B",
                                WaterbodyName == "Jennie Creek" ~ "B",
                                WaterbodyName == "Tributary to Smokey Hollow" ~ "B",
                                WaterbodyName == "Spring Brook" ~ "B",
                                WaterbodyName == "Lowry Creek" ~ "C",
                                WaterbodyName == "Lowery Creek" ~ "C",
                                WaterbodyName == "Unnamed trib to Maple Dale Creek" ~ "D",
                                WaterbodyName == "Little Pine Creek" ~ "D",
                                WaterbodyName == "Fancy Creek" ~ "D",
                                WaterbodyName == "Tamarack Creek" ~ "D",
                                WaterbodyName == "Little Scarboro Creek" ~ "D",
                                WaterbodyName == "Tagatz Creek" ~ "E",
                                WaterbodyName == "Caves Creek" ~ "E",
                                WaterbodyName == "Spranger Creek" ~ "E",
                                WaterbodyName == "Flume Creek" ~ "E",
                                WaterbodyName == "Laxey Creek" ~ "F",
                                WaterbodyName == "Melancthon Creek" ~ "F",
                                WaterbodyName == "Byrds Creek" ~ "F",
                                WaterbodyName == "Chester Creek" ~ "F",
                                WaterbodyName == "Little Silver Creek" ~ "F",
                                WaterbodyName == "Foulds Creek" ~ "F",
                                WaterbodyName == "South Fork Hay River" ~ "F",
                                WaterbodyName == "Fourmile Creek" ~ "F",
                                WaterbodyName == "Buena Vista Creek" ~ "F",
                                WaterbodyName == "Schuett Creek" ~ "F",
                                WaterbodyName == "Tomorrow River" ~ "F",
                                WaterbodyName == "Plover River" ~ "F",
                                WaterbodyName == "Tenmile Creek - South Branch" ~ "F",
                                WaterbodyName == "Devils Creek" ~ "F",
                                WaterbodyName == "Bruce Creek" ~ "F",
                                WaterbodyName == "Cady Creek" ~ "F",
                                WaterbodyName == "Walczak Creek" ~ "F",
                                WaterbodyName == "Tisch Mills Creek" ~ "F",
                                WaterbodyName == "Lunch Creek" ~ "F")) %>% 
  select(WaterbodyName, WBIC, Tree_group) %>% 
  distinct()

write_csv(Tree_groupings, 
          file = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/Trees/Tree_groups.csv")
