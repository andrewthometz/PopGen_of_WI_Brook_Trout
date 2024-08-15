library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)

# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Prep 2111 data to work with plotting
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain",
         HUC_2 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_2)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

fctLevelOrder <- Samples_2205 %>% 
  distinct(WaterbodyName) %>% 
  arrange(WaterbodyName)

Samples_2205 <- Samples_2205 %>% 
  mutate(WaterbodyName = fct_relevel(WaterbodyName, fctLevelOrder$WaterbodyName))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

###

gstMatrix <- pairwise_Gst_Nei(Data_2205)

gst_tidy <- gstMatrix %>% 
  tidy()

mean_distances <- gst_tidy %>% 
  group_by(item1) %>% 
  summarize(average_dist = round(mean(distance), digits = 3))

gst_tidy %>% write_csv("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/Gst_Fst/pwise_dist_2205.csv")

heatmap <- gst_tidy %>% 
  mutate(distance = round(distance, digits = 2)) %>% 
  ggplot(aes(x = item1, y = item2, fill = distance)) +
  geom_tile(color = "grey") +
  scale_fill_gradient(low = "blue", 
                      high = "red", 
                      space = "Lab", 
                      name = bquote("Nei's"~G[ST]),
                      breaks = c(round(min(gst_tidy$distance), 2),
                                 round(mean(gst_tidy$distance), 2),
                                 round(max(gst_tidy$distance), 2))) +
  scale_y_discrete(limits = rev(rownames(as.matrix(gstMatrix)))) +
  scale_x_discrete(position = "bottom") +
  labs(x = "",
       y = "") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -60, 
                                   vjust = 0.5, 
                                   size = 6, 
                                   hjust = 0.01,
                                   color = "black"),
        axis.text.y = element_text(size = 6,
                                   color = "black"))

ggsave(filename = "GST_heatmap.pdf",
       plot = heatmap,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 8,
       width = 10,
       units = "in")

ggsave(filename = "GST_heatmap.png",
       plot = heatmap,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 8,
       width = 10,
       units = "in")


#### Exact G test ####


# genepop package
test_diff("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_genepop.txt", 
          genic = FALSE,
          pairs = FALSE,
          outputFile = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Gst_Fst/Genepop_exact_test.txt")


#### Experimenting with other packages ####
# adegenet package
genpop_list <- as.genpop(Data_2205$tab)
dist.genpop(genpop_list)
nei.dist(Data_2205)

# hierfstat package
hierfstat_2205 <- genind2hierfstat(Data_2205)
genet.dist(hierfstat_2205)
pairwise.neifst(hierfstat_2205)
pairwise.WCfst(hierfstat_2205)

allele_freqs <- pop.freq(hierfstat_2205)
pp.fst(allele_freqs) # Couldn't get this to work