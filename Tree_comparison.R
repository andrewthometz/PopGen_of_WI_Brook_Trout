library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(radiator)
library(RColorBrewer)

library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

#### Read in Amplicon genepop file ####
Amp_gen <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Amplicon_genepop.gen", 
                        ncode = 3L, 
                        quiet = FALSE)

# Read in CE vs Amp population data
Amp_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  filter(SampleID %in% rownames(Amp_gen@tab)) %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location) %>% 
  arrange(match(SampleID, rownames(Amp_gen@tab)))

Amp_pop_data$WaterbodyName <- paste(Amp_pop_data$WaterbodyName, "_Amp", sep = "")

# Fill pop slot
Amp_gen@pop <- as_factor(Amp_pop_data$WaterbodyName)

#### Read in data ####
Data_2111 <- read.genepop("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/2111_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  arrange(Cohort) %>% 
  filter(SampleID %in% rownames(Data_2111@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2111@tab))) %>% 
  mutate(WaterbodyName = case_when(str_detect(WaterbodyName, "St. Croix") ~ "Hatchery",
                                   .default = WaterbodyName))

# Fill 2111 pop slot and subset to just domestic fish
Data_2111@pop <- as_factor(Samples_2111$WaterbodyName)

Domestics <- popsub(Data_2111, 
                    sublist = "Hatchery",
                    drop = FALSE)

# Filter to CE vs Amp loci
Domestics <- Domestics[loc = locNames(Amp_gen)]

# Read in 2205 genetic data for reference comparisons
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Filter to CE vs Amp loci
temp_gen <- Data_2205[loc = locNames(Amp_gen)]

# Write and read it back in to work around repooling error
#temp_gen %>% 
#  tidy_genind() %>% 
#  write_genepop(genepop.header = "59 Survey pops with just 7 CE vs Amp loci",
#                filename = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/2205_pops_7_loci")

Data_2205_7_loci <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/2205_pops_7_loci.gen", 
                                 ncode = 3L, 
                                 quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205_7_loci@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205_7_loci@tab))) %>% 
  select(WaterbodyName, SampleID, WBIC, HUC_8, HUC_6, HUC_4, HUC_2)

Data_2205_7_loci@pop <- as_factor(Samples_2205$WaterbodyName)

# Repool
All_dat_amplicon <- repool(Data_2205_7_loci, Domestics, Amp_gen)

#### Create tree ####
# Build initial tree (creates phylo object)
Phylo_tree <- aboot(All_dat_amplicon,
                    strata = All_dat_amplicon@pop,
                    distance = "nei.dist")

# Add HUC info to tree
temp_metadata <- Samples_2205 %>% 
  select(-SampleID) %>% 
  rename(label = WaterbodyName) %>% 
  distinct()

temp_tree <- Phylo_tree %>% 
  as_tibble() %>% 
  left_join(temp_metadata) 

# Convert phylo object into treedata object
Tree_data <- as.treedata(temp_tree)

# Plot treedata object using ggtree
tree_1 <- ggtree(Tree_data,
                 branch.length = "none") + # Use branchlength = "none" for cladogram
  geom_tiplab(size = 3.5,
              #linesize = 1,
              show.legend = FALSE) +
  scale_color_discrete(breaks = c("Chippewa",
                                  "Northwestern Lake Michigan",
                                  "Rock",
                                  "Southwestern Lake Michigan",
                                  "St. Croix",
                                  "Upper Mississippi-Maquoketa-Plum",
                                  "Western Lake Superior",
                                  "Wisconsin")) +
  theme_tree(legend.position = "bottom",
             legend.direction = "horizontal",
             plot.margin = unit(c(0, 0, 0, 0), "in"),
             legend.key.size = unit(0.25, 'in'),
             legend.key.height = unit(0.25, 'in'),
             legend.key.width = unit(0.25, 'in'),
             legend.title = element_text(size = 16), #change legend title font size
             legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title = "Watershed (HUC 4)",
                              override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Tree_amplicon.pdf",
       plot = tree_1,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 12,
       width = 30,
       units = "in")

######################## Do this again with converted CE data ##################################

#### Read in converted CE genepop file ####
CE_converted <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/CE_genepop_converted.gen", 
                             ncode = 3L, 
                             quiet = FALSE)

# Read in CE population data
CE_pop_data <- read_excel("X:/2201_BKT_msat_conversion/samplesForCEComparison.xlsx") %>% 
  filter(SampleID %in% rownames(CE_converted@tab)) %>% 
  select(SampleID, Location) %>% 
  rename(WaterbodyName = Location) %>% 
  arrange(match(SampleID, rownames(CE_converted@tab)))

CE_pop_data$WaterbodyName <- paste(CE_pop_data$WaterbodyName, "_CE", sep = "")

# Fill 2111 pop slot and subset to just domestic fish
CE_converted@pop <- as_factor(CE_pop_data$WaterbodyName)

# Repool
All_dat_converted <- repool(Data_2205_7_loci, Domestics, CE_converted)

#### Create tree ####
# Build initial tree (creates phylo object)
Phylo_tree_2 <- aboot(All_dat_converted,
                    strata = All_dat_converted@pop,
                    distance = "nei.dist")

# Add HUC info to tree
temp_metadata_2 <- Samples_2205 %>% 
  select(-SampleID) %>% 
  rename(label = WaterbodyName) %>% 
  distinct()

temp_tree_2 <- Phylo_tree_2 %>% 
  as_tibble() %>% 
  left_join(temp_metadata) 

# Convert phylo object into treedata object
Tree_data_2 <- as.treedata(temp_tree_2)

# Plot treedata object using ggtree
tree_2 <- ggtree(Tree_data_2,
                 branch.length = "none") + # Use branchlength = "none" for cladogram
  geom_tiplab(size = 3.5,
              #linesize = 1,
              show.legend = FALSE) +
  scale_color_discrete(breaks = c("Chippewa",
                                  "Northwestern Lake Michigan",
                                  "Rock",
                                  "Southwestern Lake Michigan",
                                  "St. Croix",
                                  "Upper Mississippi-Maquoketa-Plum",
                                  "Western Lake Superior",
                                  "Wisconsin")) +
  theme_tree(legend.position = "bottom",
             legend.direction = "horizontal",
             plot.margin = unit(c(0, 0, 0, 0), "in"),
             legend.key.size = unit(0.25, 'in'),
             legend.key.height = unit(0.25, 'in'),
             legend.key.width = unit(0.25, 'in'),
             legend.title = element_text(size = 16), #change legend title font size
             legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title = "Watershed (HUC 4)",
                              override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Tree_converted_CE.pdf",
       plot = tree_2,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 12,
       width = 30,
       units = "in")

######################## Make tree using just converted CE pops and Amp pops ##################################
# Repool
CE_and_Amp <- repool(Amp_gen, CE_converted)

#### Create tree ####
# Build initial tree (creates phylo object)
Phylo_tree_3 <- aboot(CE_and_Amp,
                      strata = CE_and_Amp@pop,
                      distance = "nei.dist")

# Add HUC info to tree
#temp_metadata_3 <- Samples_2205 %>% 
#  select(-SampleID) %>% 
#  rename(label = WaterbodyName) %>% 
#  distinct()

temp_tree_3 <- Phylo_tree_3 %>% 
  as_tibble() #%>% 
#  left_join(temp_metadata) 

# Convert phylo object into treedata object
Tree_data_3 <- as.treedata(temp_tree_3)

# Plot treedata object using ggtree
tree_3 <- ggtree(Tree_data_3,
                 branch.length = "none") + # Use branchlength = "none" for cladogram
  geom_tiplab(size = 3.5,
              #linesize = 1,
              show.legend = FALSE) +
  scale_color_discrete(breaks = c("Chippewa",
                                  "Northwestern Lake Michigan",
                                  "Rock",
                                  "Southwestern Lake Michigan",
                                  "St. Croix",
                                  "Upper Mississippi-Maquoketa-Plum",
                                  "Western Lake Superior",
                                  "Wisconsin")) +
  theme_tree(legend.position = "bottom",
             legend.direction = "horizontal",
             plot.margin = unit(c(0, 0, 0, 0), "in"),
             legend.key.size = unit(0.25, 'in'),
             legend.key.height = unit(0.25, 'in'),
             legend.key.width = unit(0.25, 'in'),
             legend.title = element_text(size = 16), #change legend title font size
             legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title = "Watershed (HUC 4)",
                              override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Tree_CE_and_Amp.pdf",
       plot = tree_3,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 12,
       width = 30,
       units = "in")

############################### Create tree with converted CE data, amp data, and my pops ##########################
All_data <- repool(Amp_gen, CE_converted, Domestics, Data_2205_7_loci)

#### Create tree ####
# Build initial tree (creates phylo object)
Phylo_tree_4 <- aboot(All_data,
                      strata = All_data@pop,
                      distance = "nei.dist")

# Add HUC info to tree
#temp_metadata_4 <- Samples_2205 %>% 
#  select(-SampleID) %>% 
#  rename(label = WaterbodyName) %>% 
#  distinct()

temp_tree_4 <- Phylo_tree_4 %>% 
  as_tibble() # %>% 
  #left_join(temp_metadata_4) 

# Convert phylo object into treedata object
Tree_data_4 <- as.treedata(temp_tree_4)

# Plot treedata object using ggtree
tree_4 <- ggtree(Tree_data_4,
                 branch.length = "none") + # Use branchlength = "none" for cladogram
  geom_tiplab(size = 3.5,
              #linesize = 1,
              show.legend = FALSE) +
  scale_color_discrete(breaks = c("Chippewa",
                                  "Northwestern Lake Michigan",
                                  "Rock",
                                  "Southwestern Lake Michigan",
                                  "St. Croix",
                                  "Upper Mississippi-Maquoketa-Plum",
                                  "Western Lake Superior",
                                  "Wisconsin")) +
  theme_tree(legend.position = "bottom",
             legend.direction = "horizontal",
             plot.margin = unit(c(0, 0, 0, 0), "in"),
             legend.key.size = unit(0.25, 'in'),
             legend.key.height = unit(0.25, 'in'),
             legend.key.width = unit(0.25, 'in'),
             legend.title = element_text(size = 16), #change legend title font size
             legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title = "Watershed (HUC 4)",
                              override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Tree_CE_Amp_and_2205.pdf",
       plot = tree_4,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/CE_vs_Amp/Plots",
       height = 12,
       width = 30,
       units = "in")

