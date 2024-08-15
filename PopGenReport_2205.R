library(tidyverse)
library(readxl)
#library(MCGLfxns)
library(adegenet)
library(HardyWeinberg) # HWE and LD
library(PopGenReport) # Genetic diversity stats
library(pegas) # AMOVA
library(poppr)

#### Read in genetic data ####
Data_2205 <- read.genepop("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/2205_genepop.gen", 
                          ncode = 3L, 
                          quiet = FALSE)

setwd("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/PopGenReport/")

#### Takes about 1 hour to run ####

popgenreport(
  Data_2205,
  mk.counts = TRUE,
  mk.map = FALSE,
  mk.locihz = TRUE,
  mk.hwe = TRUE,
  mk.fst = FALSE,
  mk.gd.smouse = FALSE,
  mk.gd.kosman = FALSE,
  mk.pcoa = FALSE,
  mk.spautocor = FALSE,
  mk.allele.dist = TRUE,
  mk.null.all = TRUE,
  mk.allel.rich = TRUE,
  mk.differ.stats = FALSE,
  mk.custom = FALSE,
  fname = "PGR_2205_4.23.23",
  foldername = "PGR_2205_4.23.23",
  path.pgr = getwd(),
  mk.Rcode = FALSE,
  mk.complete = FALSE,
  mk.pdf = FALSE
)
