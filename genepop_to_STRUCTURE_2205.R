library(tidyverse)
library(readxl)
library(adegenet)
library(genepop)
library(mmod)
library(hierfstat)
library(broom)
library(poppr)

#### Read genind file containing random subsample of 30 st. croix fish data ####
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/63pops_plus_30domestics.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

# Read in 2111 metadata
Samples_2111 <- read_delim("X:/2111_F1F2D_BKT/2111analysis/Thometz_scripts/Samples_2111.csv") %>% 
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab))) %>% 
  mutate(WaterbodyName = case_when(str_detect(WaterbodyName, "St. Croix") ~ "St. Croix Falls domestic",
                                   .default = WaterbodyName)) %>% 
  select(SampleID, WaterbodyName)

# Read in 2205 metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv") %>% 
  bind_rows(Samples_2111) %>%
  filter(SampleID %in% rownames(Data_2205@tab)) %>% 
  arrange(match(SampleID, rownames(Data_2205@tab))) 

Data_2205@pop <- as_factor(Samples_2205$WaterbodyName)

#### Convert data to STRUCTURE #### (These did not work well, use function below instead)
library(graph4lg)

#test <- popsub(Data_2205, sublist = c("Tomorrow River", "Bruce Creek"))

#genind_to_structure(test, output = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/STRUCTURE/test.txt")

#genind_to_structure(Joined_2205_2111, output = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/STRUCTURE/STR_input_2205_2111.txt")

#### Convert data to STRUCTURE ####
devtools::install_github("rystanley/genepopedit") 
library(genepopedit)

#### Trying this function ####
genind2structure <- function(obj, file = "", pops = FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind = rep(adegenet::indNames(obj), each = pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each = pl)
    tab <- cbind(tab, data.frame(pop = popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow = dim(tab)[1], ncol = adegenet::nLoc(obj),
                           dimnames = list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == adegenet::indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
}

#genind2structure(Data_2205, 
#                 pops = FALSE,
#                 file = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/STRUCTURE_input_2205.txt")

genind2structure(Data_2205,
                 pops = FALSE,
                 file = "X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Structure_relatedness/STRUCTURE/Final_run/STRUCTURE_final_input_2205.txt")
