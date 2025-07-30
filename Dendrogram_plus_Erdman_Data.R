# Load packages
library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
#library(RColorBrewer)
library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

#######################################################################################################################
#### Build trees to visualize genetic relatedness among my survey populations and Brad Erdman's survey populations ####
#######################################################################################################################

# Read in Master brook trout genepop file
UNIFIED_BKT <- read.genepop("X:/filepath.../Erdman_integration/UNIFIED_BKT_genepop.gen",
                            ncode = 3L,
                            quiet = FALSE)

# Read in 2111 project metadata
Samples_2111 <- read_delim("X:/filepath.../Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls domestic") %>% 
  select(SampleID, WaterbodyName)

# Read in 2205 project metadata
Samples_2205 <- read_delim("X:/filepath.../Samples_2205.csv") %>% 
  bind_rows(Samples_2111)

# Read in Brad Erdman's genotype data
Erdman_samples <- read_excel("X:/filepath.../Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx")

# Bind the metadata
All_metadata <- Samples_2205 %>% 
  bind_rows(Erdman_samples) %>% 
  select(SampleID, WBIC, WaterbodyName, HUC_4, HUC_8) %>% 
  filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) 

# Assign pop slot
UNIFIED_BKT@pop <- as_factor(All_metadata$WaterbodyName)

#### Build the tree ####
# Build the initial tree (creates a phylo object)
Phylo_tree <- aboot(UNIFIED_BKT,
                    strata = UNIFIED_BKT@pop,
                    distance = "nei.dist",
                    tree = "nj") # Do "nj" instead of default "upgma" to make dendrogram

# Get HUC info to color-code leaf tips
temp_metadata <- All_metadata %>% 
  rename(label = WaterbodyName) %>% 
  select(-SampleID) %>% 
  distinct()

# Turn the phylo object into a tibble to add HUC data
tree_tibble <- Phylo_tree %>% 
  as_tibble() %>% 
  left_join(temp_metadata)

# Convert the tibble into a treedata object for ggtree plotting
Tree_data <- as.treedata(tree_tibble)

# Plot the treedata object using ggtree (dendrogram)
tree_1 <- ggtree(Tree_data, 
                 aes(color = HUC_4), 
                 size = 1) +
  geom_tiplab(show.legend = FALSE,
              size = 2.5) +
  geom_treescale(color = "grey40",
                 linesize = 0.5) +
  scale_color_discrete(name = "Subregion (HUC 4)",
                       na.value = "grey40",
                       breaks = c("Chippewa",
                                  "Northwestern Lake Michigan",
                                  "Rock",
                                  "Southwestern Lake Michigan",
                                  "St. Croix",
                                  "Upper Mississippi-Black-Root",
                                  "Upper Mississippi-Maquoketa-Plum",
                                  "Western Lake Superior",
                                  "Wisconsin")) +
  xlim(0, 0.75) + # This can help make tree fit
  theme_tree(legend.position = c(0.8, 0.8),
             legend.key.size = unit(0.25, 'in'),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Tree_UNIFIED_BKT.png",
       plot = tree_1,
       device = "png",
       path = "X:/filepath.../Erdman_integration/Plots_figures",
       height = 12,
       width = 9,
       units = "in")
