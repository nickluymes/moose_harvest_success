######################################################
############ Moose Tag Fill Rates Project ############ 
######################################################
### Understanding what factors impact hunter success rates and developing
### a model to predict success rates for upcoming hunting seasons

######################################################
### Script to combine and prepare data for use in Bayesian hierarchical 
### survival models 

library(tidyverse)
source("set_up.R")

# Filter data for Time to Event analysis ----

tidy_tte <- readRDS("data/processed/tidy_long.rds") %>%
  left_join(readRDS("data/processed/tidy_season.rds")) %>%
  ### keep only hunters that had a lottery tag, hunted in their tag's WMU 
  ### for at least 1 day and less than the total number of days in the season
  filter(
    Tag_Lottery == TRUE,
    Hunted == TRUE,
    Tag_Match_WMU == TRUE,
    Days > 0,
    Days <= Season_Length,
    !is.na(Tag_FA),
    !is.na(Tag_Moose)
  ) %>%
  mutate(
    ### successful hunters where those that harvested a moose in a WMU that 
    ### matched the WMU and age/sex class of their tag
    ### unsuccessful hunters did not harvest a moose with their tag
    Tag_Success = as.integer(
      case_when(
        UsedTag & Tag_Match_Moose & Success_Match_WMU ~ "1",
        UsedTag == FALSE & is.na(Tag_Match_Moose) & 
          is.na(Success_Match_WMU) ~ "0",
        TRUE ~ NA_character_
      )
    ),
    ### binned days for success model and ending hunt model
    ### binned days for ending hunt model every multiple needs to be higher
    ### than the highest binned day for success for the conditional probs 
    ### to work
   # Days2 = ifelse(Tag_Success==1, ceiling(Days/5), floor(Days/5))
  ) %>%
  filter(!is.na(Tag_Success), Tag_Lottery,
         ### filter out bins that span less than 5 days 
        # Days2*5<=Season_Length, Days2>0
        )

tidy_tte %>%
  ggplot(aes(x=Days)) +
  geom_histogram(binwidth=1) +
  facet_wrap(~Tag_Success)

### create dataset for if hunters ended their hunt without being successful
tidy_stop <- tidy_tte %>%
  ### success = stopping hunt without successful harvest
  mutate(Tag_stop = ifelse(Tag_Success == 0,1,0))

### create dataset of harvest success during each day a hunter hunted
### 0 = no success, 1 = success
tidy_tte_1 <- tidy_tte %>%
  as.data.frame() %>%
  ### this package automates the process of adding extra rows for the number of 
  ### days until success of censorship (stopping hunt)
  discSurv::dataLong(
    timeColumn = "Days",
    eventColumn = "Tag_Success",
    timeAsFactor = FALSE
  ) %>%
  as_tibble() %>%
  ### calculate the number of successful hunters during each day of a hunt
  group_by(Year, Tag_WMU, Tag_Res, Tag_FA, Tag_Moose, 
           timeInt, Season_Start, Season_Length) %>%
  summarise(y = sum(y), n = n())

### create dataset of hunters ending hunts during each day a hunter hunted
### 0 = hunted next day, 1 = ended hunt
tidy_stop_tte <- tidy_stop %>%
  as.data.frame() %>%
  discSurv::dataLong(
    timeColumn = "Days",
    eventColumn = "Tag_stop",
    timeAsFactor = FALSE
  ) %>%
  as_tibble() %>%
  ### calculate the number of hunters that finished their hunt without a moose
  ### during each day of a hunt
  group_by(Year, Tag_WMU, Tag_Res, Tag_FA, Tag_Moose, 
           timeInt, Season_Start, Season_Length) %>%
  summarise(y = sum(y), n = n()) %>%
  ### change the total number of hunters each day to the total number of hunters minus
  ### the number of hunters that were successful because they do not contribute to 
  ### the pool of hunters that could end their hunt before the next day *without*
  ### harvesting a moose
  left_join(tidy_tte_1 %>% dplyr::select(Year, Tag_WMU, Tag_Res, Tag_FA, Tag_Moose, 
                                         timeInt, Season_Start, Season_Length, y) %>%
              rename(ysuc = y)) %>%
  mutate(n = n - ysuc,
         ysuc = NULL) %>%
  ### keep only those days with some unsuccessful hunters remaining
  filter(n > 0)



# Combine Datasets ----
### combine covariate and time to event data and transform covariates into
### spatial and temporal components for time to successful harvest data
tidy_success <- readRDS("data/processed/tidy_tags.rds") %>%
  full_join(readRDS("data/processed/tidy_precip_daily.rds")) %>%
  full_join(readRDS("data/processed/tidy_temp_daily.rds")) %>%
  full_join(readRDS("data/processed/tidy_snow_weeks.rds") %>%
              select(Year, Tag_WMU, Tag_Moose, Tag_Res,
                            Tag_FA, snow)) %>%
  full_join(readRDS("data/processed/tidy_mai.rds") %>%
              select(Year, Tag_WMU, Tag_Moose, cpop_dens, cpop_size_harv)) %>%
  full_join(readRDS("data/processed/tidy_hunters.rds")) %>%
  full_join(readRDS("data/processed/tidy_road.rds") %>%
              dplyr::select(Tag_WMU, road_density)) %>%
  ### calculate WMU-means and deviations from means
  group_by(Tag_WMU, Tag_Moose, Tag_FA, Tag_Res) %>%
  mutate(pop_avg = mean(cpop_dens, na.rm = TRUE),
         pop_change = cpop_dens - pop_avg,
         dhunters_avg = mean(dhunters, na.rm = TRUE),
         dhunters_change = dhunters - dhunters_avg,
         ### divide tags issued by # of moose and log transform
         tagpermooselog = log(Tag_Issued/cpop_size_harv),
         tagpermooselog_avg = mean(tagpermooselog, na.rm = TRUE),
         tagpermooselog_change = tagpermooselog - tagpermooselog_avg,
         precip_avg = mean(precip, na.rm = TRUE),
         precip_change = precip - precip_avg,
         temp_avg = mean(temp, na.rm = TRUE),
         temp_change = temp - temp_avg,
         snow_avg = mean(snow, na.rm = TRUE),
         snow_change = snow - snow_avg) %>%
  ### join hunter time to event data
  right_join(tidy_tte_1) %>%
  mutate(Year = as.factor(Year)) %>%
  ### keep only entries for adult moose with mai and tags issued data
  filter(Tag_Moose != "Calf", 
         !is.na(cpop_dens),
         !is.na(Tag_Issued)) %>%
  ungroup() %>%
  mutate(
    ### ensure categories match grouping levels
    Tag_FA = fct_relevel(Tag_FA, c("G","A")),
    Tag_WMU = fct_relevel(Tag_WMU, WMU_level),
    ### create polynomial transformation for day of hunt
    timeInt_poly = poly(timeInt, 3),
    timeInt_poly1 = timeInt_poly[,1],
    timeInt_poly2 = timeInt_poly[,2],
    timeInt_poly3 = timeInt_poly[,3],
    ### trying log transforming timeint before taking polynomials
    timeInt_log = log(timeInt),
    timeInt_log_poly = poly(timeInt_log, 3),
    timeInt_log_poly1 = timeInt_log_poly[,1],
    timeInt_log_poly2 = timeInt_log_poly[,2],
    timeInt_log_poly3 = timeInt_log_poly[,3],
    ### set road density for 7A to zero - no major roads identified in that area
    road_density = ifelse(is.na(road_density), 0, road_density)
  )
saveRDS(tidy_success, "data/processed/tidy_success.rds")

### combine covariate and time to event data and transform covariates into
### spatial and temporal components for time to ending hunt data
tidy_unsuccess <- readRDS("data/processed/tidy_tags.rds") %>%
  full_join(readRDS("data/processed/tidy_precip_daily.rds")) %>%
  full_join(readRDS("data/processed/tidy_temp_daily.rds")) %>%
  full_join(readRDS("data/processed/tidy_snow_weeks.rds") %>%
              select(Year, Tag_WMU, Tag_Moose, Tag_Res,
                     Tag_FA, snow)) %>%
  full_join(readRDS("data/processed/tidy_mai.rds") %>%
              select(Year, Tag_WMU, Tag_Moose, cpop_dens, cpop_size)) %>%
  full_join(readRDS("data/processed/tidy_hunters.rds")) %>%
  full_join(readRDS("data/processed/tidy_road.rds") %>%
              dplyr::select(Tag_WMU, road_density)) %>%
  ### calculate WMU-means and deviations from means
  group_by(Tag_WMU, Tag_Moose, Tag_FA, Tag_Res) %>%
  mutate(pop_avg = mean(cpop_dens, na.rm = TRUE),
         pop_change = cpop_dens - pop_avg,
         dhunters_avg = mean(dhunters, na.rm = TRUE),
         dhunters_change = dhunters - dhunters_avg,
         ### divide tags issued by # of moose and log transform
         tagpermooselog = log(Tag_Issued/cpop_size),
         tagpermooselog_avg = mean(tagpermooselog, na.rm = TRUE),
         tagpermooselog_change = tagpermooselog - tagpermooselog_avg,
         precip_avg = mean(precip, na.rm = TRUE),
         precip_change = precip - precip_avg,
         temp_avg = mean(temp, na.rm = TRUE),
         temp_change = temp - temp_avg,
         snow_avg = mean(snow, na.rm = TRUE),
         snow_change = snow - snow_avg) %>%
  ### join hunter time to event data
  right_join(tidy_stop_tte) %>%
  mutate(Year = as.factor(Year)) %>%
  ### keep only entries for adult moose with mai and tags issued data
  filter(Tag_Moose != "Calf", 
         !is.na(cpop_dens),
         !is.na(Tag_Issued)) %>%
  ungroup() %>%
  mutate(
    ### ensure categories match grouping levels
    Tag_FA = fct_relevel(Tag_FA, c("G","A")),
    Tag_WMU = fct_relevel(Tag_WMU, WMU_level),
    ### create polynomial transformation for day of hunt
    timeInt_poly = poly(timeInt, 3),
    timeInt_poly1 = timeInt_poly[,1],
    timeInt_poly2 = timeInt_poly[,2],
    timeInt_poly3 = timeInt_poly[,3],
    ### trying log transforming timeint before taking polynomials
    timeInt_log = log(timeInt),
    timeInt_log_poly = poly(timeInt_log, 3),
    timeInt_log_poly1 = timeInt_log_poly[,1],
    timeInt_log_poly2 = timeInt_log_poly[,2],
    timeInt_log_poly3 = timeInt_log_poly[,3],
    ### set road density for 7A to zero - no major roads identified in that area
    road_density = ifelse(is.na(road_density), 0, road_density)
  )
saveRDS(tidy_unsuccess, "data/processed/tidy_unsuccess.rds")

tidy_unsuccess %>%
  ggplot(aes(x=y)) +
  geom_histogram(binwidth=1) +
  xlim(-2,20)

# Create adjacency matrix ----
### use spdep package to create a matrix with 1s and 0s where a 1 indicates 
### that the WMU in the corresponding row and column are neighbours.
adjacency_mat <- readRDS("data/processed/WMU_shp.rds") %>%
  spdep::poly2nb() %>%
  spdep::nb2mat(style = "B")
rownames(adjacency_mat) <- WMU_shp$WMU
saveRDS(adjacency_mat, "data/processed/adjacency_mat.rds")

# tidy_cez <- tibble(Tag_WMU = WMU_level_order) %>%
#   mutate(
#     wmu_num = as.numeric(gsub(".*?([0-9]+).*", "\\1", Tag_WMU)),
#     CEZ = case_when(
#       wmu_num %in% c(1,2,16,17,18,24,25,26,27) ~ "A",
#       wmu_num %in% c(3,4,14,15,19,21,22,23,30,33,34) ~ "B",
#       wmu_num %in% c(5,8,12,13) ~ "C1",
#       Tag_WMU %in% c("9A","11B") ~ "C1",
#       wmu_num %in% c(28,29,31,32,35,38,39,40,41) ~ "C2",
#       wmu_num %in% c(6,7) ~ "D1",
#       Tag_WMU %in% c("9B","11A") ~ "D1",
#       wmu_num %in% c(36,37,42,46,47,48,49,50,53,54,55,56,57,58,60,61,62,63) ~ "D2",
#       wmu_num %in% c(59,65) ~ "E3"
#     ),
#     wmu_num = NULL
#   )

