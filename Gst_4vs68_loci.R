library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)

# Run GST at just 4 loci to see if they match results from 68 loci
#### Read in data ####
# Read in 2205 genetic data at just 4 loci of interest
Data_2205_4loci <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/2205_63pops_4loci.gen", 
                                ncode = 3L, 
                                quiet = FALSE)

# 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls domestic",
         HUC_2 = "Hatchery",
         HUC_4 = "Hatchery",
         HUC_8 = "Hatchery",
         .keep = "unused") %>% 
  select(SampleID, WaterbodyName, HUC_8, HUC_4, HUC_2)

# 2205 metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>% 
  filter(SampleID %in% rownames(Data_2205_4loci@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205_4loci@tab)))

fctLevelOrder <- Samples_2205 %>% 
  distinct(WaterbodyName) %>% 
  arrange(WaterbodyName)

Samples_2205 <- Samples_2205 %>% 
  mutate(WaterbodyName = fct_relevel(WaterbodyName, fctLevelOrder$WaterbodyName))

# Fill the pop slots
Data_2205_4loci@pop <- as_factor(Samples_2205$WaterbodyName)

# Run the analysis
gstMatrix <- pairwise_Gst_Nei(Data_2205_4loci)

gst_tidy <- gstMatrix %>% 
  tidy()

gst_tidy %>% write_csv("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Locus_testing/Gst_4loci.csv")

# Plot on heatmap
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

ggsave(filename = "GST_heatmap_4loci.pdf",
       plot = heatmap,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 10,
       units = "in")

ggsave(filename = "GST_heatmap_4loci.png",
       plot = heatmap,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 10,
       units = "in")

#### Join gst estimates with 68 loci estimates to compare #### This is implemented in the "Diversity_4vs68_loci.R" script

gst_68loci <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Gst_Fst/pwise_Nei_Gst_plus2111.csv") %>% 
  filter(item1 != "St. Croix domestic",
         item2 != "St. Croix domestic")

gst_comparison <- gst_tidy %>% 
  rename(dist_4loci = distance) %>% 
  left_join(gst_68loci) %>% 
  rename(Pop_1 = item1,
         Pop_2 = item2,
         dist_68loci = distance) %>% 
  mutate(dist_4loci = round(dist_4loci, 3),
         dist_68loci = round(dist_68loci, 3))

# Run linear regression to compare
lm_dist = lm(dist_4loci ~ dist_68loci, data = gst_comparison)

summary(lm_dist)

gst_comparison %>% 
  mutate(dist_4loci = case_when(dist_4loci < 0 ~ 0,
                                .default = dist_4loci)) %>% 
  ggplot(aes(x = dist_68loci, y = dist_4loci)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  annotate(geom = "text", 
           label = "(1:1)",
           x = 0.2, 
           y = 0.22,
           color = "red") +
  annotate(geom = "text",
           label = bquote(R^2~"= 0.499"),
           x = 0.22, 
           y = 0.14) +
  scale_y_continuous(limits = c(0, 0.25),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 0.25),
                     expand = c(0,0)) +
  labs(x = "68 loci",
       y = "4 loci",
       title = bquote("Nei's"~G[ST])) +
  theme_classic()

  


