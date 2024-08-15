library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(radiator)

# Run genetic diversity at just 4 loci to see if they match results from 68 loci

#### Read in data ####
# Read in 2205 genetic data at just 4 loci of interest
Data_2205_4loci <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/2205_63pops_4loci.gen", 
                                ncode = 3L, 
                                quiet = FALSE)

# 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(Cohort == "Domestic") %>% 
  mutate(WaterbodyName = "St. Croix Falls Strain",
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

# Fill the pop slots
Data_2205_4loci@pop <- as_factor(Samples_2205$WaterbodyName)

#### Calculate genetic diversity measures ####
library(hierfstat)

gd <- basic.stats(Data_2205_4loci)

H_expected <- as.data.frame.matrix(gd$Hs) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:63,
               names_to = "pop",
               values_to = "He") %>% 
  drop_na(He) %>% 
  group_by(pop) %>% 
  summarize("He" = round(mean(He), 3))

H_observed <- as.data.frame.matrix(gd$Ho) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:63,
               names_to = "group",
               values_to = "Ho") %>% 
  drop_na(Ho) %>% 
  group_by(group) %>% 
  summarize("Ho" = round(mean(Ho), 3))

F_is <- as.data.frame.matrix(gd$Fis) %>% 
  as_tibble(rownames = "locus") %>% 
  pivot_longer(2:63,
               names_to = "group",
               values_to = "Fis") %>% 
  drop_na(Fis) %>% 
  group_by(group) %>% 
  summarize("Fis" = round(mean(Fis), 3))

Allelic_richness <- as.data.frame(allelic.richness(Data_2205_4loci)) %>% 
  as_tibble(rownames = "locus") %>% 
  select(-min.all) %>% 
  pivot_longer(2:63,
               names_to = "group",
               values_to = "ar") %>% 
  drop_na(ar) %>% 
  group_by(group) %>% 
  summarize("ar" = round(mean(ar), 3))

GD_4loci <- bind_cols(H_expected$pop,
                      Allelic_richness$ar, 
                      H_observed$Ho, 
                      H_expected$He,
                      F_is$Fis) %>% 
  rename("WaterbodyName" = 1,
         "Ar_4loci" = 2,
         "Ho_4loci" = 3, 
         "He_4loci" = 4,
         "Fis_4loci" = 5)

# Read in gd data from 68 loci
GD_68loci <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Genetic_diversity/Diversity_2205.csv") %>% 
  #select(-c(1, WBIC, Ne)) %>% 
  rename(Ar_68loci = Ar,
         Ho_68loci = Ho,
         He_68loci = He,
         Fis_68loci = Fis)

# Join to compare
GD_comparison <- GD_4loci %>% 
  left_join(GD_68loci)

# Plot Ar
lm_Ar <- lm(Ar_4loci ~ Ar_68loci, data = GD_comparison)

summary(lm_Ar)

Ar_plot <- GD_comparison %>% 
  ggplot(aes(x = Ar_68loci, y = Ar_4loci)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  # annotate(geom = "text",
  #          label = "(1:1)",
  #          x = 1.35,
  #          y = 4,
  #          color = "red") +
  annotate(geom = "text",
           label = bquote({R^{2}}[adj]~"= 0.525"),
           x = 1.34, 
           y = 4.2) +
  annotate(geom = "text",
           label = bquote(p-value~"= 1.694"~{"e"^{-11}}),
           x = 1.34, 
           y = 4) +
  scale_y_continuous(limits = c(NA, 5),
                     expand = c(0,0.001)) +
  scale_x_continuous(limits = c(NA, 1.5),
                     expand = c(0,0.001)) +
  labs(x = "",
       y = "",
       title = bquote("(B)  Allelic richness"~(A[r]))) +
  theme_classic() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Plot Ho
lm_Ho <- lm(Ho_4loci ~ Ho_68loci, data = GD_comparison)

summary(lm_Ho)

Ho_plot <- GD_comparison %>% 
  ggplot(aes(x = Ho_68loci, y = Ho_4loci)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  annotate(geom = "text",
           label = "(1:1)",
           x = 0.45,
           y = 0.48,
           color = "red") +
  annotate(geom = "text",
           label = bquote({R^{2}}[adj]~"= 0.190"),
           x = 0.33, 
           y = 0.65) +
  annotate(geom = "text",
           label = bquote(p-value~"= 0.0002"),
           x = 0.33, 
           y = 0.62) +
  scale_y_continuous(#limits = c(0, NA),
                     expand = c(0,0.001)) +
  scale_x_continuous(#limits = c(0, NA),
                     expand = c(0,0.001)) +
  labs(x = "68 loci",
       y = "4 loci",
       title = bquote("(C)  Observed heterozygosity"~(H[o]))) +
  theme_classic() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Plot He
lm_He <- lm(He_4loci ~ He_68loci, data = GD_comparison)

summary(lm_He)

He_plot <- GD_comparison %>% 
  ggplot(aes(x = He_68loci, y = He_4loci)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  annotate(geom = "text",
           label = "(1:1)",
           x = 0.48,
           y = 0.5,
           color = "red") +
  annotate(geom = "text",
           label = bquote({R^{2}}[adj]~"= 0.346"),
           x = 0.33, 
           y = 0.68) +
  annotate(geom = "text",
           label = bquote(p-value~"= 2.953"~{"e"^{-7}}),
           x = 0.33, 
           y = 0.66) +
  scale_y_continuous(#limits = c(0, NA),
                     expand = c(0,0.001)) +
  scale_x_continuous(#limits = c(0, NA),
                     expand = c(0,0.001)) +
  labs(x = "68 loci",
       y = "",
       title = bquote("(D)  Expected heterozygosity"~(H[e]))) +
  theme_classic() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Plot Fis (Did not use Fis for final)
lm_Fis <- lm(Fis_4loci ~ Fis_68loci, data = GD_comparison)

summary(lm_Fis)

Fis_plot <- GD_comparison %>% 
  ggplot(aes(x = Fis_68loci, y = Fis_4loci)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm", color = "blue") +
  annotate(geom = "text",
           label = "(1:1)",
           x = 0.09,
           y = 0.06,
           color = "red") +
  annotate(geom = "text",
           label = bquote(R^2~"= 0.140"),
           x = 0.01, 
           y = 0.1) +
  scale_y_continuous(#limits = c(0, NA),
    expand = c(0,0)) +
  scale_x_continuous(#limits = c(0, NA),
    expand = c(0,0)) +
  labs(x = "68 loci",
       y = "4 loci",
       title = bquote("Inbreeding coefficient"~(F[IS]))) +
  theme_classic()

#### Bring in Gst data for four panel plot ####
gst_4loci <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Locus_testing/Gst_4loci.csv")

gst_68loci <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/Gst_Fst/pwise_dist_2205.csv")

gst_comparison <- gst_4loci %>% 
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

gst_plot <- gst_comparison %>% 
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
           label = bquote({R^{2}}[adj]~"= 0.494"),
           x = 0.2, 
           y = 0.04) +
  annotate(geom = "text",
           label = bquote(p-value~"= 2.2"~{"e"^{-16}}),
           x = 0.2, 
           y = 0.023) +
  scale_y_continuous(limits = c(0, 0.25),
                     expand = c(0,0.001)) +
  scale_x_continuous(limits = c(0, 0.25),
                     expand = c(0,0.001)) +
  labs(x = "",
       y = "4 loci",
       title = bquote("(A)  Pairwise genetic distance")) +
  theme_classic() +
  theme(plot.margin = margin(10, 10, 10, 10))

#### Now arrange them nicely in four panel plot ####
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

comparison_plot <- (gst_plot | Ar_plot) / (Ho_plot | He_plot)

ggsave(filename = "Loci_4vs68_plot.pdf",
       plot = comparison_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 10,
       units = "in")

ggsave(filename = "Loci_4vs68_plot.png",
       plot = comparison_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Erdman_integration/Plots_figures",
       height = 8,
       width = 10,
       units = "in")
