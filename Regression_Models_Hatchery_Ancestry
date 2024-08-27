# Load necessary packages
library(car)
library(ggsignif)
library(tidyverse)
library(readxl)
library(adegenet)
#library(betareg)
#library(gamlss) # This package allows for zero-inflated beta regression (BEZI)
#library(zoib)
library(brms)

##########################################################################################
#### Run regression models to identify potential relationships with hatchery ancestry ####
##########################################################################################

# Read in project metadata
Samples_2205 <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Samples_2205.csv")

# Read in WDNR stocking data for relevant HUC12s
Stocking_data <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Stocking_histories_2205.csv")

# Read in genetic diversity data (contains Ne estimates)
Diversity <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Genetic_diversity/Diversity_2205.csv") %>% 
  select(-1)

# Read in hatchery introgression data
Hatchery_ID <- read_delim("X:/2205_BKT_feral_broodstock_ID/Thometz_scripts/Analyses/Hatchery_introgression/Hatchery_ID_2205.csv") %>% 
  select(-1)

# Create stocking variables at HUC12 level
HUC12_stocking <- Stocking_data %>% 
  group_by(HUC_12) %>% 
  summarize(Total_n_stocked = sum(n_stocked),
            n_stocking_events = n(),
            Yrs_since_stocking = round(mean(Stocking_Year), 0))

# Join data together to create df for modeling
Modeling_data <- Samples_2205 %>% 
  select(WBIC, WaterbodyName, HUC_12) %>% 
  distinct() %>% 
  left_join(HUC12_stocking) %>% 
  left_join(Diversity) %>% 
  #left_join(Ne_estimates) %>% 
  #mutate(Ne = as.numeric(Ne)) %>% 
  left_join(Hatchery_ID) %>% 
  mutate(Total_n_stocked = replace_na(Total_n_stocked, 0)) %>% 
  mutate(n_stocking_events = replace_na(n_stocking_events, 0)) 

##################################
#### Describe model variables ####
##################################

#### Response variables: ####
# Proportion of fish with hatchery assignment probability >= 25% (scrapped this because not very useful)
# Proportion of fish with hatchery assignment probability >= 75% (scrapped this becuase not very useful)
# Average assignment probability to hatchery ancestry

#### Predictor variables: ####
# Total number of fish stocked in HUC12 watershed
# Number of stocking events in HUC12 watershed
# Years since the mean of stocking years in HUC12 watershed
# Ne among stocked populations
# Relative stocking intensity (e.g., number of fish stocked / Ne) ?

#### Average assignment probability to hatchery ancestry ####
# Hypothesis 1:
# Total number of fish stocked in HUC12 watershed is positively related to mean hatchery identity
Modeling_data_temp <- Modeling_data %>% 
  select(-Yrs_since_stocking) # Doing this because NAs in this column are problematic

model_1a <- brm(Ave_hatchery_ID ~ Total_n_stocked,
                family = zero_inflated_beta(),
                data = Modeling_data)

summary(model_1a)

Modeling_data %>% 
  ggplot(aes(x = Total_n_stocked, y = Ave_hatchery_ID)) + 
  geom_smooth(method = "glm") +
  geom_point() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  xlab("Total number of fish stocked in HUC 12 watershed") +
  ylab("Mean hatchery identity") +
  theme_classic()

# Hypothesis 2: 
# Number of stocking events in HUC12 watershed is positively related to mean hatchery identity
model_1b <- gamlss(Ave_hatchery_ID ~ n_stocking_events, 
                   family = BEZI, 
                   data = na.omit(Modeling_data))

summary(model_1b)

Modeling_data %>% 
  ggplot(aes(x = n_stocking_events, y = Ave_hatchery_ID)) + 
  geom_smooth(method = "glm") +
  geom_point() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  xlab("Number of stocking events in HUC 12 watershed") +
  ylab("Mean hatchery identity") +
  theme_classic()

# Hypothesis 3: 
# Mean stocking year in HUC12 watershed is positively related to mean hatchery identity
model_1c <- gamlss(Ave_hatchery_ID ~ Yrs_since_stocking, 
                   family = BEZI, 
                   data = na.omit(Modeling_data))

summary(model_1c)

Modeling_data %>% 
  ggplot(aes(x = Yrs_since_stocking, y = Ave_hatchery_ID)) + 
  geom_smooth(method = "glm") +
  geom_point() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  xlab("Mean stocking year") +
  ylab("Mean hatchery identity") +
  theme_classic()

# Hypothesis 4: 
# Ne is negatively related to mean hatchery identity
Modeling_data2 <- Modeling_data %>% 
  filter(Ne != "Inf")

model_1d <- gamlss(Ave_hatchery_ID ~ Ne, 
                family = BEZI, 
                data = na.omit(Modeling_data2))

summary(model_1d)

Modeling_data2 %>% 
  ggplot(aes(x = Ne, y = Ave_hatchery_ID)) + 
  geom_smooth(method = "glm") +
  geom_point() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  xlab("Effective population size (Ne)") +
  ylab("Mean hatchery identity") +
  theme_classic()
