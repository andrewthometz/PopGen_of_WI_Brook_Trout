library(tidyverse)
library(ggforce)

### Function to make length_dist files tidy ###

length_dist_tidy <- function(x){
  
  BKT_ID <- deparse(substitute(x))
  
  output <- x %>% 
    mutate(uSat_locus = str_replace_all(Microsatellite, "-", "_"), .keep = "unused") %>% 
    select(-sum) %>%
    pivot_longer(-c(uSat_locus, scores), names_to = "Length", values_to = "Read_count") %>% 
    add_column(SampleID = BKT_ID) %>% 
    mutate(SampleID = str_remove_all(SampleID, "BKT_")) %>% 
    mutate(Length = as.numeric(Length))
  
  print(output)
}

### Read in 12 BKT ### Grabbed largest few files from each cat folder. Used some from each folder to get "even" representation

BKT_18_06622 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_1of5_Msat_output/length_distribution/Genotype_18-06622.txt")
BKT_18_06180 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_1of5_Msat_output/length_distribution/Genotype_18-06180.txt")
BKT_15_03167 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_1of5_Msat_output/length_distribution/Genotype_15-03167.txt")

BKT_21_10773 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_2of5_Msat_output/length_distribution/Genotype_21-10773.txt")
BKT_19_21562 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_2of5_Msat_output/length_distribution/Genotype_19-21562.txt")
BKT_18_07013 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_2of5_Msat_output/length_distribution/Genotype_18-07013.txt")

BKT_22_12362 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_3of5_Msat_output/length_distribution/Genotype_22-12362.txt")
BKT_22_10464 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_3of5_Msat_output/length_distribution/Genotype_22-10464.txt")

BKT_22_12613 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_4of5_Msat_output/length_distribution/Genotype_22-12613.txt")
BKT_22_12710 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_4of5_Msat_output/length_distribution/Genotype_22-12710.txt")

BKT_22_13858 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_5of5_Msat_output/length_distribution/Genotype_22-13858.txt")
BKT_22_15856 <- read_delim("X:/2205_BKT_feral_broodstock_ID/2205_MEGAsat_outputs/2205_paired_5of5_Msat_output/length_distribution/Genotype_22-15856.txt")

### Run each fish's length_dist file through my tidying function while joining them ###

All_LD_data <- bind_rows(length_dist_tidy(BKT_18_06622),
                         length_dist_tidy(BKT_18_06180),
                         length_dist_tidy(BKT_15_03167),
                         length_dist_tidy(BKT_21_10773),
                         length_dist_tidy(BKT_19_21562),
                         length_dist_tidy(BKT_18_07013),
                         length_dist_tidy(BKT_22_12362),
                         length_dist_tidy(BKT_22_10464),
                         length_dist_tidy(BKT_22_12613),
                         length_dist_tidy(BKT_22_12710),
                         length_dist_tidy(BKT_22_13858),
                         length_dist_tidy(BKT_22_15856))

All_LD_data %>% 
  #filter(Read_count > 0) %>% 
  group_by(uSat_locus, SampleID) %>% 
  #slice_max(Read_count, n = 10) %>% 
  filter(Read_count >= max(Read_count)*0.25) %>% 
  ungroup() %>% 
  summarize(x = n_distinct(uSat_locus))

### Make peak morphology plots, one locus and 12 BKT per page ### Saves directly as pdf

# Beware, this takes about 40min to run #

pdf("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_PeakMorph_Thometz.pdf", paper = "a4r", width = 11, height = 9)

ProgressBar <- txtProgressBar(min = 0, max = 91, style = 3)

for(i in 1:91){
  
print(All_LD_data %>% 
        group_by(uSat_locus, SampleID) %>% 
        #slice_max(Read_count, n = 8) %>% ### This seems to determine run time (4 ~ 40min, 10 ~ 1.5hrs)
        filter(Read_count >= max(Read_count)*0.05) %>% 
        ggplot(aes(x = Length, y = Read_count)) +
        geom_col(fill = "seagreen") +
        geom_text(aes(label = Length), 
                  hjust = 0.5, 
                  vjust = 1.5, 
                  colour = "black") +
        ylab("Read count") +
        xlab("Allele length") +
        scale_y_continuous(expand = c(0,0)) +
        #scale_x_continuous(breaks = c(0:200), expand = c(0,0)) +
        theme_classic() +
        facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(SampleID), as.factor(scores)), 
                            nrow = 3,
                            ncol = 4,
                            page = i, 
                            scales = "free"))
  Sys.sleep(0.1)
  setTxtProgressBar(ProgressBar, i)
}
close(ProgressBar)
dev.off()

### Single page of peak morphs for testing plot changes ###

All_LD_data %>% 
  #filter(Read_count > 0) %>% 
  group_by(uSat_locus, SampleID) %>% 
  #slice_max(Read_count, n = 10) %>% 
  filter(Read_count >= max(Read_count)*0.05) %>% 
  ggplot(aes(x = Length, y = Read_count)) +
  geom_col(fill = "seagreen") +
  geom_text(aes(label = Length), hjust = 0.5, vjust = 1.5, colour = "black") +
  ylab("Read count") +
  xlab("Allele length") +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(breaks = c(0:200), expand = c(0,0)) +
  theme_classic() +
  facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(SampleID), as.factor(scores)), 
                      nrow = 3,
                      ncol = 4,
                      page = 49, 
                      scales = "free")

#### Produce plot of peak morphs for 9 loci, 1 fish #### (For suppelemental figures)
loci_list <- All_LD_data %>% 
  distinct(uSat_locus)

Filtered_tibble <- All_LD_data %>% 
  group_by(uSat_locus) %>% 
  filter(SampleID == "22_10464",
         Read_count >= max(Read_count)*0.05,
         uSat_locus == "SFO_18" |
         uSat_locus == "L_SFOC28" |
         uSat_locus == "SFOC86" |
         uSat_locus == "SfoD91" |
         uSat_locus == "L_Ssa_13.6" |
         uSat_locus == "L_Ssa_16.2" |
         uSat_locus == "Salv_2_tri_20_30" |
         uSat_locus == "Sfon_6_Di_01_402" |
         uSat_locus == "Sfon_6_tri_20_20") %>% 
  mutate(uSat_locus = case_when(uSat_locus == "SFO_18" ~ "Sfo-18",
                                uSat_locus == "L_SFOC28" ~ "L-SFOC28",
                                uSat_locus == "SFOC86" ~ "SfoC86",
                                #uSat_locus == "" ~ "",
                                uSat_locus == "L_Ssa_13.6" ~ "L-Ssa-13.6",
                                uSat_locus == "L_Ssa_16.2" ~ "L-Ssa-16.2",
                                uSat_locus == "Salv_2_tri_20_30" ~ "Salv-2_tri-20-30",
                                uSat_locus == "Sfon_6_Di_01_402" ~ "Sfon-6_Di-01-402",
                                uSat_locus == "Sfon_6_tri_20_20" ~ "Sfon-6_tri-20-20",
                                .default = uSat_locus))

# Function to only display integers axes
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

# Plot it
onefish_nineloci <- Filtered_tibble %>% 
  ggplot(aes(x = Length, y = Read_count)) +
  geom_col(fill = "grey55") +
  geom_text(aes(label = Length), hjust = 0.5, vjust = 1.5, colour = "black") +
  ylab("Read count") +
  xlab("Allele length") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = integer_breaks()) +
  theme_classic() +
  facet_wrap_paginate(facets = vars(as.factor(uSat_locus), as.factor(scores)), 
                      nrow = 3,
                      ncol = 3,
                      page = 1, 
                      scales = "free")

ggsave(filename = "Peakmorph_1fish_9loci.png",
       plot = onefish_nineloci,
       device = "png",
       path = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Polished_plots_figures",
       height = 8,
       width = 8,
       units = "in")




