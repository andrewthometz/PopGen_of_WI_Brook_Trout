library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
#library(RColorBrewer)
library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

#### Read in Master brook trout genepop file ####
UNIFIED_BKT <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/UNIFIED_BKT_genepop.gen",
                            ncode = 3L,
                            quiet = FALSE)

#### Read in metadata ####
# 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls domestic") %>% 
  select(SampleID, WaterbodyName)

# 2205 metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111)

# Erdman
Erdman_samples <- read_excel("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Erdman_WI_BKT_Genotypes.xlsx")

# Bind the metadata
All_metadata <- Samples_2205 %>% 
  bind_rows(Erdman_samples) %>% 
  select(SampleID, WBIC, WaterbodyName, HUC_4, HUC_8) %>% 
  filter(SampleID %in% rownames(UNIFIED_BKT@tab)) %>% 
  arrange(match(SampleID, rownames(UNIFIED_BKT@tab))) 

# Assign pop slot
UNIFIED_BKT@pop <- as_factor(All_metadata$WaterbodyName)

#### Build the tree ####
# Build initial tree (creates phylo object)
Phylo_tree <- aboot(UNIFIED_BKT,
                    strata = UNIFIED_BKT@pop,
                    distance = "nei.dist",
                    tree = "nj") # Do "nj" instead of default "upgma" to make dendrogram

# Get HUC info to color code leaf tips
temp_metadata <- All_metadata %>% 
  rename(label = WaterbodyName) %>% 
  select(-SampleID) %>% 
  distinct()

# Turn phylo object into tibble to add huc data
tree_tibble <- Phylo_tree %>% 
  as_tibble() %>% 
  left_join(temp_metadata)

# Convert tibble into treedata object for ggtree plotting
Tree_data <- as.treedata(tree_tibble)

# Plot treedata object using ggtree (dendrogram)
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
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 12,
       width = 9,
       units = "in")

###########################################



# Plot treedata object using ggtree
tree_1 <- ggtree(Tree_data, aes(color = HUC_4),
                 branch.length = "none",
                 layout = "circular") + # Use branchlength = "none" for cladogram
  geom_tiplab(size = 3.5,
              #linesize = 1,
              show.legend = FALSE) +
  # scale_color_discrete(breaks = c("Chippewa",
  #                                 "Northwestern Lake Michigan",
  #                                 "Rock",
  #                                 "Southwestern Lake Michigan",
  #                                 "St. Croix",
  #                                 "Upper Mississippi-Maquoketa-Plum",
  #                                 "Western Lake Superior",
  #                                 "Wisconsin")) +
  #scale_color_manual(values=c("red", "blue", "green", "orange", "yellow", "pink", "brown", "purple")) +
  #xlim(-1,1) +
  #geom_treescale(x = 0, y = 45) +
  #geom_tippoint() +
  #geom_nodepoint(color = "black") +
  #geom_hilight() +
  #geom_range() +
  #geom_cladelab() +
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

ggsave(filename = "Tree_ALL_BKT.pdf",
       plot = tree_1,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 16,
       width = 16,
       units = "in")
