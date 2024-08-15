#### Run Mantel correlogram to evaluate relationship between geo dist and gen dist ####

########################################################
######         Create Mantel correlograms         ######
########################################################

##### Load required packages #####
library(vegan)
library(mmod)
#install_github("jaredhomola/VPLandscapeGenetics")
library(VPLandscapeGenetics)
library(tidyverse)
library(geosphere)
data(VPLandscapeGenetics)

# Read in 2205 genetic data
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_All_63pops.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab)))

# Fill the pop slots
Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

##### Calculate genetic distances #####
nei.distances <- pairwise_Gst_Nei(Data_2205, linearized = TRUE)
nei.dist.matrix <- as.matrix(nei.distances)
dimnames(nei.dist.matrix) <- list(1:63, 1:63)

##### Calculate geographic distances #####
Coords <- Samples_2205 %>% 
  select(Longitude, Latitude) %>%
  distinct() 

geo.dist.matrix <- distm(Coords, fun = distGeo) %>% as.matrix()
geo.distances <- geo.dist.matrix/1000
dim(geo.distances) <- c(63, 63)
dimnames(geo.distances) <- list(1:63, 1:63)

###### Run test and produce default plot ######
break.pts <- c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360)

correlogram <- mantel.correlog(nei.dist.matrix, 
                               D.geo = geo.distances, 
                               mult = "BH", 
                               nperm = 1000, 
                               break.pts = break.pts, 
                               cutoff = TRUE)

# Basic plot
plot(correlogram, alpha = 0.05)

## Assign pch based on corrected p value using a new column
bkt.results <- correlogram$mantel.res %>% 
  as.data.frame() %>% 
  drop_na() %>% 
  rename("p.corrected" = "Pr(corrected)") %>% 
  mutate(Significance = case_when(p.corrected <= 0.05 ~ "Significant",
                                  .default = "Non-significant"))

##### Publication plot #####
pub.plot <- bkt.results %>% 
  ggplot(aes(x = class.index, y = Mantel.cor)) +
  geom_hline(yintercept = 0, 
             color = "red") +
  geom_line() +
  geom_point(aes(fill = as_factor(Significance)),
             color = "black",
             shape = 21, # shape = 21 is imperative for color and fill to work together
             show.legend = FALSE,
             size = 3) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.03, 0.08)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 200)) +
  labs(x = "Distance class (km)",
       y = "Mantel correlation coefficient") +
  theme_classic() +
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(filename = "Mantel_correlation.pdf",
       plot = pub.plot,
       device = "pdf",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 4,
       width = 5,
       units = "in")

ggsave(filename = "Mantel_correlation.png",
       plot = pub.plot,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 4,
       width = 5,
       units = "in")
