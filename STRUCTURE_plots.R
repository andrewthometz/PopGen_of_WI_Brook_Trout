library(tidyverse)
library(readxl)
library(adegenet)
library(poppr)
library(ggh4x)
library(patchwork)

# Read in 2205 sample data
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  select(WaterbodyName, SampleID, WBIC, HUC_2, HUC_4, HUC_6, HUC_8)

# K = 2 STRUCTURE run
K2 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/STRUCTURE/Reruns/K=2/K2_Results/K2_Results.txt") %>% 
  select(-n) %>% 
  mutate(C2 = as.numeric(C2)) %>% 
  pivot_longer(cols = c(C1, C2), 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205)

C1 <- K2 %>% 
  group_by(SampleID) %>% 
  filter(Probability == max(Probability)) %>% 
  ungroup() %>% 
  group_by(WBIC) %>% 
  count(Cluster) %>% 
  filter(n == max(n),
         Cluster == "C1") %>% 
  select(-c(Cluster, n)) %>% 
  mutate(WBIC = as_factor(WBIC))

C2 <- K2 %>% 
  group_by(SampleID) %>% 
  filter(Probability == max(Probability)) %>% 
  ungroup() %>% 
  group_by(WBIC) %>% 
  count(Cluster) %>% 
  filter(n == max(n),
         Cluster == "C2") %>% 
  select(-c(Cluster, n)) %>% 
  mutate(WBIC = as_factor(WBIC))

K2 <- K2 %>% 
  mutate(pop_clust = case_when(WBIC %in% C1$WBIC ~ "C1",
                               WBIC %in% C2$WBIC ~ "C2"))

# Initial K=2 structure plot for hierarchical structure figure
K2_plot <- K2 %>% 
  ggplot(aes(SampleID, Probability, fill = factor(Cluster))) +
  geom_col() +
  facet_grid(cols = vars(pop_clust, HUC_8, WaterbodyName), 
             switch = "x", 
             scales = "free", 
             space = "free") +
  labs(x = "", y = "Assignment probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.25, 'in'),
        legend.key.height = unit(0.25, 'in'),
        legend.key.width = unit(0.25, 'in'),
        legend.title = element_text(size = 24), # change legend title font size
        legend.text = element_text(size = 20)) + 
  guides(fill = guide_legend(title = "Genetic lineage"))

ggsave(filename = "Initial_K2_str_plot.pdf",
       plot = K2_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Hierarchical_structure_plots",
       height = 4,
       width = 20,
       units = "in")

# K = 11 STRUCTURE run
K11 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/STRUCTURE/Reruns/K=11/K11_Results/K11_Results.txt") %>% 
  select(-1) %>% 
  mutate(C11 = as.numeric(C11)) %>% 
  pivot_longer(cols = C1:C11, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205)

K11 %>% 
  ggplot(aes(SampleID, Probability, fill = factor(Cluster))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~HUC_8 + WaterbodyName, 
             switch = "x", 
             scales = "free", 
             space = "free") +
  labs(x = "Waterbody", y = "Assignment probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_discrete("Dark2") +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank()) + 
  guides(fill = guide_legend(title = "Genetic lineage")) +
  theme(strip.text.x = element_text(angle = -90))

# K = 13 STRUCTURE run
K13 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Reruns/K=13/K13_Results/K13_Results.txt") %>% 
  select(-1) %>% 
  mutate(C13 = as.numeric(C13)) %>% 
  pivot_longer(cols = C1:C13, 
               names_to = "Cluster", 
               values_to = "Probability") %>% 
  left_join(Samples_2205)

K13_plot <- K13 %>% 
  ggplot(aes(SampleID, Probability, fill = factor(Cluster))) +
  geom_col(color = "grey", linewidth = 0.1) +
  facet_grid(cols = vars(HUC_8), 
             switch = "x", 
             scales = "free", 
             space = "free") +
  labs(x = "HUC_8", 
       y = "Assignment probability",
       fill = "Genetic cluster") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #scale_fill_discrete("Dark2") +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        #strip.text.x = element_text(angle = -90)) +
        strip.text.x = element_text()) 

ggsave(filename = "K13_plot.pdf",
       plot = K13_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 10,
       width = 10,
       units = "in")

ggsave(filename = "K13_plot.png",
       plot = K13_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Presentation_plots",
       height = 2,
       width = 30,
       units = "in")

# Zoomed in
zoomed <- K13 %>% 
  filter(HUC_8 == "Brule" |
         HUC_8 == "Lake Dubay")

K13_zoomed_plot <- zoomed %>% 
  ggplot(aes(SampleID, Probability, fill = factor(Cluster))) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(cols = vars(HUC_8), 
             switch = "x", 
             scales = "free", 
             space = "free") +
  labs(x = "HUC_8", 
       y = "Assignment probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  #scale_fill_discrete("Dark2") +
  theme_minimal(base_size = 15) + 
  theme(panel.spacing.x = unit(0.01, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        #strip.text.x = element_text(angle = -90)) +
        strip.text.x = element_text()) + 
  guides(fill = guide_legend(title = "Genetic cluster"))

ggsave(filename = "K13_zoomed_plot.png",
       plot = K13_zoomed_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures/Presentation_plots",
       height = 3,
       width = 7,
       units = "in")

