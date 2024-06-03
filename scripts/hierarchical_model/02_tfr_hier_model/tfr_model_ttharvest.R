######################################################
############ Moose Tag Fill Rates Project ############ 
######################################################
### Understanding what factors impact hunter success rates and developing
### a model to predict success rates for upcoming hunting seasons

######################################################
### Script to run Bayesian Hierarchical Survival Model for time to successful
### hunt data. This script is intended to be run on compute canada servers

library(tidyverse)
library(brms)

tidy_success <- readRDS("data/processed/tidy_success.rds") %>%
  group_by(Tag_WMU, Tag_Res, Tag_FA, Tag_Moose, timeInt) %>%
  arrange(Year) %>%
  mutate(Season_Start_change = Season_Start - mean(Season_Start))

temp <- readRDS("data/processed/adjacency_mat.rds")
### keep only WMUs in the tidy_success data set
mat_ids <- which(rownames(temp)%in%tidy_success$Tag_WMU)
adj_mat <- temp[mat_ids,mat_ids]

### set priors for beta coefficients
b_prior <- set_prior("normal(0, 1)", class = "b")
mod <- brm(y | trials(n) ~ scale(timeInt_poly1) +
            scale(timeInt_poly2) +
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
             (1 + Tag_Res*Tag_Moose*Tag_FA|Tag_WMU) + 
             (1|Year) +
              car(M, gr = Tag_WMU, type = "icar"),
           family = binomial(link = "cloglog"),
           data = tidy_success, 
           prior = b_prior,
           save_pars = save_pars(all = TRUE),
           data2 = list(M = adj_mat), chains = 3, 
           iter = 4000, cores = 3,
           seed = 5, init = "0",
           control = list(adapt_delta = 0.8),
           sample_prior = "only")
saveRDS(mod, "brms_ttharvest_mod.rds")
