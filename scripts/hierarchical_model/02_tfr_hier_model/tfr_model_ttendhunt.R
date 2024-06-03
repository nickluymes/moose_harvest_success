######################################################
############ Moose Tag Fill Rates Project ############ 
######################################################
### Understanding what factors impact hunter success rates and developing
### a model to predict success rates for upcoming hunting seasons

######################################################
### Script to run Bayesian Hierarchical Survival Model for time to ending
### hunt data. This script is intended to be run on compute canada servers

library(tidyverse)
library(brms)

tidy_unsuccess <- readRDS("tidy_unsuccess.rds") %>%
  filter(timeInt < Season_Length) %>%
  ### calculate probability of hunt day lining up with the end of the season
  mutate(pdays = timeInt/Season_Length)

temp <- readRDS("adjacency_mat.rds")
### keep only WMUs in the tidy_success data set
mat_ids <- which(rownames(temp)%in%tidy_unsuccess$Tag_WMU)
adj_mat <- temp[mat_ids,mat_ids]


b_prior <- set_prior("normal(0, 1)", class = "b")
mod <- brm(y | trials(n) ~ scale(timeInt_log_poly1) +
             scale(timeInt_log_poly2) +
             scale(timeInt_log_poly3) +
             scale(pdays) +
             scale(tagpermooselog_avg) +
             scale(tagpermooselog_change) +
             scale(pop_avg) +
             scale(pop_change) +
             scale(Season_Length) +
             scale(precip_change) +
             scale(temp_change) +
             scale(snow_avg) +
             scale(snow_change) +
             scale(road_density) +
             scale(dhunters_avg) +
             scale(dhunters_change) +
             Tag_Res*Tag_Moose*Tag_FA +
             (1|Tag_WMU) + 
             (1|Year),
           family = binomial(link = "cloglog"),
           data = tidy_unsuccess, 
           prior = b_prior,
           data2 = list(M = adj_mat), chains = 4, 
           iter = 4000, cores = 4,
           seed = 5, inits = "0",
           control = list(adapt_delta = 0.8))
saveRDS(mod, "brms_mod_ttendhunt.rds")

