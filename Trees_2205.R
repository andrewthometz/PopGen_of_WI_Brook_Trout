library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(RColorBrewer)

library(BiocManager)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

#### Build trees! ####
# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Prep 2111 data to work with plotting
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain") %>% 
  select(SampleID, WaterbodyName)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

# Build initial tree (creates phylo object)
Phylo_tree <- aboot(Data_2205,
                    strata = Data_2205@pop,
                    distance = "nei.dist",
                    tree = "nj",
                    cutoff = 50) # Do "nj" instead of default "upgma" to make dendrogram

# Get HUC info to color code leaf tips
temp_metadata <- Samples_2205 %>% 
  rename(label = WaterbodyName) %>% 
  select(-SampleID) %>% 
  distinct()

# Turn phylo object into tibble to add huc data
tree_tibble <- Phylo_tree %>% 
  as_tibble() %>% 
  left_join(temp_metadata) %>% 
  mutate(#Cluster = case_when(label %in% Samples_revised$Cluster ~ label, .default = NA),
         Bootstraps = case_when(!(label %in% Samples_2205$WaterbodyName) ~ label, .default = NA))

# Convert tibble into treedata object for ggtree plotting
Tree_data <- as.treedata(tree_tibble)

# Plot treedata object using ggtree (dendrogram)
brewer.pal(n = 8, name = "Set1")

my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow3", "#A65628", "#F781BF")

my_breaks <- temp_metadata %>% 
  select(HUC_4) %>% 
  distinct() %>% 
  drop_na()

tree_1 <- ggtree(Tree_data, 
                 aes(color = HUC_4), 
                 size = 1) +
  geom_tiplab(show.legend = FALSE) +
  geom_treescale(color = "black",
                 linesize = 1) +
  geom_text(aes(label = Bootstraps),
            color = "black",
            hjust = 1.2,
            vjust = -0.5,
            size = 3,
            show.legend = FALSE) +
  geom_strip(taxa1 = "Venison Creek", 
             taxa2 = "Noisy Creek",
             label = "A",
             align = FALSE,
             barsize = 1,
             hjust = -0.4,
             offset = 0.128,
             extend = 0.35,
             fontsize = 4) +
  geom_strip(taxa1 = "Lepage Creek", 
             taxa2 = "Spring Brook",
             label = "B",
             align = FALSE,
             barsize = 1,
             hjust = -0.4,
             offset = 0.086,
             extend = 0.35,
             fontsize = 4) +
  geom_strip(taxa1 = "Lowry Creek", 
             taxa2 = "Lowery Creek",
             label = "C",
             align = FALSE,
             barsize = 1,
             hjust = -0.4,
             offset = 0.03,
             extend = 0.35,
             fontsize = 4) +
  geom_strip(taxa1 = "Unnamed trib to Maple Dale Creek", 
             taxa2 = "Little Scarboro Creek",
             label = "D",
             align = FALSE,
             barsize = 1,
             hjust = -0.4,
             offset = 0.076,
             extend = 0.35,
             fontsize = 4) +
  geom_strip(taxa1 = "Tagatz Creek", 
             taxa2 = "Flume Creek",
             label = "E",
             align = FALSE,
             barsize = 1,
             hjust = -0.4,
             offset = 0.036,
             extend = 0.35,
             fontsize = 4) +
  geom_strip(taxa1 = "Laxey Creek", 
             taxa2 = "Lunch Creek",
             label = "F",
             align = FALSE,
             barsize = 1,
             hjust = -0.4,
             offset = 0.065,
             extend = 0.35,
             fontsize = 4) +
  scale_colour_manual(name = "Subregion (HUC 4)", # Can use this or scale_color_discrete()
                      na.value = "black",
                      values = my_colors,
                      breaks = my_breaks$HUC_4) +
  xlim(0, 0.25) + # This can help make tree fit
  theme_tree(legend.position = c(0.85, 0.58),
             legend.key.size = unit(0.25, 'in'),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 11)) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, linetype = 1))) 

ggsave(filename = "Best_NJtree_2205.pdf",
       plot = tree_1,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 10,
       width = 9,
       units = "in")

ggsave(filename = "Best_NJtree_2205.png",
       plot = tree_1,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Trees",
       height = 10,
       width = 9,
       units = "in")

