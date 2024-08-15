#### Modeling ####
library(car)
library(ggsignif)
library(tidyverse)
library(readxl)
library(adegenet)

# Read in master gd file
GD_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Genetic_diversity/Diversity_2205.csv") %>% 
  select(-1) #%>% 
  #filter(Study == "Thometz",
  #       WaterbodyName != "St.Croix domestic")

# Ho
Ho_mod_lat <- lm(Ho ~ Latitude, data = GD_2205)
par(mfrow = c(2,2))
plot(Ho_mod_lat)
summary(Ho_mod_lat)

# He
He_mod_lat <- lm(He ~ Latitude, data = GD_2205)
par(mfrow = c(2,2))
plot(He_mod_lat)
summary(He_mod_lat)

# Ar
Ar_mod_lat <- lm(Ar ~ Latitude, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_lat)
summary(Ar_mod_lat)

library(broom)
tidy(Ar_mod_lat)

Ar_lat_plot <- GD_2205 %>% 
  ggplot(aes(x = Latitude, y = Ar)) + 
  geom_smooth(method = "glm") +
  geom_point() +
  xlab("Latitude") +
  ylab("Allelic richness") +
  scale_x_continuous(expand = c(0,0.02)) +
  scale_y_continuous(expand = c(0,0.002)) +
  annotate(geom = "text",
           label = bquote(R^2~"= 0.191"),
           x = 45, 
           y = 1.36) +
  annotate(geom = "text",
           label = bquote(p-value~"= 0.0002"),
           x = 45, 
           y = 1.35) +
  theme_classic()

ggsave(filename = "Ar_latitude_regression.pdf",
       plot = Ar_lat_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 4.5,
       width = 5,
       units = "in")

ggsave(filename = "Ar_latitude_regression.png",
       plot = Ar_lat_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 4.5,
       width = 5,
       units = "in")

# Try the same for Ne now
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  select(SampleID, WaterbodyName, Latitude)

Ne_tidy <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Genetic_diversity/Ne/Ne_All_63pops_LD_CleanedUp.txt") %>% 
  select(SampleID, Ne) %>% 
  left_join(Samples_2205) %>% 
  select(-SampleID) %>% 
  filter(Ne > 0)

Ne_mod_lat <- lm(Ne ~ Latitude, data = Ne_tidy)
par(mfrow = c(2,2))
plot(Ne_mod_lat)
summary(Ne_mod_lat)

Ne_lat_plot <- Ne_tidy %>% 
  filter(Ne < 1000) %>% 
  ggplot(aes(x = Latitude, y = Ne)) + 
  geom_smooth(method = "lm") +
  geom_point() +
  xlab("Latitude") +
  ylab(bquote(N[e])) +
  scale_x_continuous(expand = c(0,0.02)) +
  scale_y_continuous(expand = c(0,0.2)) +
  annotate(geom = "text",
           label = bquote(R^2~"= 0.059"),
           x = 44, 
           y = 600) +
  annotate(geom = "text",
           label = bquote(p-value~"= 0.033"),
           x = 44, 
           y = 550) +
  theme_classic()

ggsave(filename = "Ne_latitude_regression.pdf",
       plot = Ne_lat_plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 4.5,
       width = 5,
       units = "in")

ggsave(filename = "Ne_latitude_regression.png",
       plot = Ne_lat_plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 4.5,
       width = 5,
       units = "in")


########################## Not using this stuff #####################################################################
# HUC 2
Ar_mod_HUC2 <- glm(Ar ~ HUC_2, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_HUC2)
summary(Ar_mod_HUC2)

GD_2205 %>% 
  ggplot(aes(x = HUC_2, y = Ar)) +
  geom_boxplot() +
  xlab("HUC 2") +
  ylab("Allelic richness") +
  geom_signif(comparisons = list(c("Great Lakes Region", "Upper Mississippi Region")), map_signif_level = TRUE) +
  theme_classic()

# HUC 4
Ar_mod_HUC4 <- glm(Ar ~ HUC_4, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_HUC4)
summary(Ar_mod_HUC4)

GD_2205 %>% 
  ggplot(aes(x = HUC_4, y = Ar)) +
  geom_boxplot() +
  xlab("HUC 4") +
  ylab("Allelic richness") +
  theme_classic()

# HUC 6
Ar_mod_HUC6 <- glm(Ar ~ HUC_6, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_HUC6)
summary(Ar_mod_HUC6)

GD_2205 %>% 
  ggplot(aes(x = HUC_6, y = Ar)) +
  geom_boxplot() +
  xlab("HUC 6") +
  ylab("Allelic richness") +
  theme_classic()

# HUC 8
Ar_mod_HUC8 <- glm(Ar ~ HUC_8, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_HUC8)
summary(Ar_mod_HUC8)

GD_2205 %>% 
  ggplot(aes(x = HUC_8, y = Ar)) +
  geom_boxplot() +
  xlab("HUC 8") +
  ylab("Allelic richness") +
  theme_classic()

#### Stocking variables ####
library(betareg)

# Total number of fish stocked
Ar_mod_nstocked <- glm(Ar ~ Total_n_stocked, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_nstocked)
summary(Ar_mod_nstocked)

GD_2205 %>% 
  ggplot(aes(x = Total_n_stocked, y = Ar)) +
  geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Total number of fish stocked") +
  ylab("Allelic richness") +
  theme_classic()
#
Fis_mod_nstocked <- glm(Fis ~ Total_n_stocked, data = GD_2205)
par(mfrow = c(2,2))
plot(Fis_mod_nstocked)
summary(Fis_mod_nstocked)

GD_2205 %>% 
  ggplot(aes(x = Total_n_stocked, y = Fis)) +
  geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Total number of fish stocked") +
  ylab("Inbreeding coefficient") +
  theme_classic()
#
Ho_mod_nstocked <- betareg(Ho ~ Total_n_stocked, data = GD_2205)
par(mfrow = c(2,2))
plot(Ho_mod_nstocked)
summary(Ho_mod_nstocked)

GD_2205 %>% 
  ggplot(aes(x = Total_n_stocked, y = Ho)) +
  #geom_smooth(method = "glm", method.args = list(family = "beta")) +
  geom_point() +  
  xlab("Total number of fish stocked") +
  ylab("Observed heterozygosity") +
  theme_classic()

# Number of stocking events
Ar_mod_stockevents <- glm(Ar ~ n_stocking_events, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_stockevents)
summary(Ar_mod_stockedevents)

GD_2205 %>% 
  ggplot(aes(x = n_stocking_events, y = Ar)) +
  geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Number of stocking events") +
  ylab("Allelic richness") +
  theme_classic()
#
Fis_mod_stockevents <- glm(Fis ~ n_stocking_events, data = GD_2205)
par(mfrow = c(2,2))
plot(Fis_mod_stockevents)
summary(Fis_mod_stockevents)

GD_2205 %>% 
  ggplot(aes(x = n_stocking_events, y = Fis)) +
  geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Number of stocking events") +
  ylab("Inbreeding coefficient") +
  theme_classic()
#
Ho_mod_stockevents <- betareg(Ho ~ n_stocking_events, data = GD_2205)
par(mfrow = c(2,2))
plot(Ho_mod_stockevents)
summary(Ho_mod_stockevents)

GD_2205 %>% 
  ggplot(aes(x = n_stocking_events, y = Ho)) +
  #geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Number of stocking events") +
  ylab("Observed heterozygosity") +
  theme_classic()

# Number of years since mean of stocking event years
Ar_mod_year <- glm(Ar ~ Years_since_stocking, data = GD_2205)
par(mfrow = c(2,2))
plot(Ar_mod_year)
summary(Ar_mod_year)

GD_2205 %>% 
  ggplot(aes(x = Years_since_stocking, y = Ar)) +
  geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Years since mean stocking year") +
  ylab("Allelic richness") +
  theme_classic()
#
Fis_mod_year <- glm(Fis ~ Years_since_stocking, data = GD_2205)
par(mfrow = c(2,2))
plot(Fis_mod_year)
summary(Fis_mod_year)

GD_2205 %>% 
  ggplot(aes(x = Years_since_stocking, y = Fis)) +
  geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Years since mean stocking year") +
  ylab("Inbreeding coefficient") +
  theme_classic()
#
Ho_mod_year <- betareg(Ho ~ Years_since_stocking, data = GD_2205)
par(mfrow = c(2,2))
plot(Ho_mod_year)
summary(Ho_mod_year)

GD_2205 %>% 
  ggplot(aes(x = Years_since_stocking, y = Ho)) +
  #geom_smooth(method = "glm") +
  geom_point() +  
  xlab("Years since mean stocking year") +
  ylab("Observed heterozygosity") +
  theme_classic()

#### Trying with multiple variable in one model #### 
Ho_mod_all <- betareg(Ho ~ Years_since_stocking + Total_n_stocked + n_stocking_events, data = GD_2205)
vif(Ho_mod_all)
par(mfrow = c(2,2))
plot(Ho_mod_year)
summary(Ho_mod_all)
