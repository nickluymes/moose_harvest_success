######################################################
############ Moose Tag Fill Rates Project ############ 
######################################################
### Understanding what factors impact hunter success rates and developing
### a model to predict success rates for upcoming hunting seasons

######################################################
### Script to analyse results from Bayesian hierarchical 
### survival models 

source("set_up.R")
library(tidyverse)
library(posterior)
library(brms)
library(tidybayes)

### load models for time to harvest and time to end hunt
mod1 <- readRDS("output/hierarchical_model/brms_modtagpermoose.rds")
mod2 <- readRDS("output/hierarchical_model/brms_unsuccess.rds")

### load associated data sets
tidy_success <- readRDS("data/processed/tidy_success.rds")
tidy_unsuccess <- readRDS("data/processed/tidy_unsuccess.rds")

### Posterior Predictive Checks

tidy_unsuccess %>%
  mutate(Tag_Moose = droplevels(Tag_Moose),
         Tag_WMU = droplevels(Tag_WMU),
         Year = droplevels(Year),
         pdays=timeInt/Season_Length) %>%
  add_predicted_draws(mod2, ndraws=1000) %>%
  group_by(.draw, .prediction, Tag_WMU) %>%
  summarise(n_rep = n()) %>%
  ungroup() %>%
  complete(.draw, .prediction=0:10, Tag_WMU, fill = list(n_rep=0)) %>%
  ggplot(aes(x = .prediction)) +
  geom_histogram(data = tidy_unsuccess, aes(x = y), binwidth = 1, fill = "lightblue") +
  stat_pointinterval(aes(y = n_rep),.width = c(.50, .80, .95, .99), size =0.3) +
  facet_wrap(~Tag_WMU, scales="free_y") +
  xlim(-2,10) +
  xlab("Number of hunts ended") +
  ylab("Frequency") +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj

ggsave("ppcheck_endinghunt.jpg", height = 8, width = 11, bg = "white")


tidy_success %>%
  mutate(Tag_Moose = droplevels(Tag_Moose),
         Tag_WMU = droplevels(Tag_WMU),
         Year = droplevels(Year),
         pdays=timeInt/Season_Length) %>%
  add_predicted_draws(mod1, ndraws=1000) %>%
  group_by(.draw, .prediction, Tag_WMU) %>%
  summarise(n_rep = n()) %>%
  ungroup() %>%
  complete(.draw, .prediction=0:10, Tag_WMU, fill = list(n_rep=0)) %>%
  ggplot(aes(x = .prediction)) +
  geom_histogram(data = tidy_success, aes(x = y), binwidth = 1, fill = "lightblue") +
  stat_pointinterval(aes(y = n_rep),.width = c(.50, .80, .95, .99), size =0.3) +
  facet_wrap(~Tag_WMU, scales="free_y") +
  xlim(-2,10) +
  xlab("Number of successful hunts") +
  ylab("Frequency") +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj

ggsave("ppcheck_huntsuccess.jpg", height = 8, width = 11, bg = "white")

### Tjur's R squared based on leave-one-out residuals

loo_r21 <- loo_R2(mod1)
loo_r22 <- loo_R2(mod2)
saveRDS(loo_r21, "output/hierarchical_model_looR2_mod1.rds")
saveRDS(loo_r22, "output/hierarchical_model_looR2_mod2.rds")


# Continuous Effects ----

### function for recursively multiplying probabilities across hunting days
### used for calculating the probability of success/ending hunt for each day 
### of the hunt
rvar_mult <- function(times, preds){
  accum = vector(mode = "list", length = 1)
  for(i in 1:length(times)) {
    if(i == 1){acc = 1 - preds[i]}
    else{
      acc = (1 - preds[i]) %**% (1 - accum[[i-1]])
    }
    accum[[i]] = 1 - acc
  }
  return(accum)
}


### function for recursively multiplying probabilities across hunting days
### alternative method for calculating the probability of success/ending hunt
### for each day of the hunt
rvar_mult2 <- function(times, preds){
  accum = vector(mode = "list", length = 1)
  vec2 = vector(mode = "list", length = 1)
  for(i in 1:length(times)) {
    if(i == 1){acc = preds[i]; vec = 1 - acc}
    else{
      acc = preds[i] %**% (vec2[[i-1]])
      vec = (1 - preds[i]) %**% (vec2[[i-1]])
    }
    accum[[i]] = acc
    vec2[[i]] = vec
  }
  return(accum)
}

rvar_func <- function(times, preds, pc){
  y_preds = 0
  new_n = 1
  for(i in 1:length(times)) {
    num = new_n
    y_pred = preds[i] %**% num
    new_n = (1 - pc[i]) %**% (num - y_pred)
    y_preds = y_preds + y_pred
  }
  return(y_preds)
}

rvar_func2 <- function(times, n, preds, pc){
  rvar_binom <- rfun(rbinom)
  y_preds = 0
  new_n = n
  for(i in 1:length(times)) {
    num = new_n
    y_pred = rvar_binom(1,num,preds[i])
    new_n = rvar_binom(1,num - y_pred,1 - pc[i])
    y_preds = y_preds + y_pred
  }
  return(y_preds)
}

cdf_effects2 <- function(param, mod, mod2, xlab, logscale = FALSE) {
  range <- modelr::seq_range(tidy_success %>% dplyr::select(sym(param)), 20)
  p1 <- tidy_success %>%
    group_by(Tag_Moose, Tag_Res, Tag_FA) %>%
    select(-timeInt, -timeInt_poly, -timeInt_poly1, -timeInt_poly2, -n, -y) %>%
    distinct() %>%
    modelr::data_grid(
      pop_avg = mean(pop_avg),
      tagpermooselog_avg = mean(tagpermooselog_avg),
      pop_change = 0,
      tagpermooselog_change = 0,
      Season_Length = 40,
      snow_avg = mean(snow_avg),
      precip_change = 0,
      temp_change = 0,
      snow_change = 0,
      road_density = mean(road_density, na.rm = TRUE),
      dhunters_avg = mean(dhunters_avg),
      dhunters_change = 0,
      timeInt = 1:Season_Length,
      Year = "1",
      n = 1, 
      Tag_WMU = "13"
    ) %>%
    filter(Tag_Res == "Res",
           Tag_Moose == "Bull",
           Tag_FA == "G") %>%
    mutate(
      Tag_Moose = as.character(Tag_Moose),
      timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
      timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
      timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3],
      pdays = timeInt/Season_Length
    ) %>%
    select(-sym(param)) %>%
    expand_grid(!!sym(param) := range) %>%
    tidybayes::add_epred_rvars(mod, allow_new_levels = TRUE, re_formula = NA) %>%
    rename(.epred1 = .epred) %>%
    tidybayes::add_epred_rvars(mod2, allow_new_levels = TRUE, re_formula = NA) %>%
    rename(.epred2 = .epred) %>%
    mutate(pop_change = 100 * pop_change,
           pop_avg = 100 * pop_avg,
           road_density = road_density/1000) %>%
    group_by(Tag_WMU, Tag_Moose, Tag_Res, Tag_FA, !!sym(param)) %>%
    arrange(timeInt) %>%
    summarise(
      cdf = rvar_func(timeInt, .epred1, .epred2),
    ) %>%
    rename(x = sym(param)) %>%
    group_by(Tag_Moose, Tag_Res, Tag_FA, x) %>%
    summarise(
      cdf = rvar_mean(cdf)
    ) %>%
    mutate(x = ifelse(rep(logscale,n()),exp(x),x)) %>%
    ggplot(aes(x = x)) +
    ggdist::stat_dist_lineribbon(aes(dist = cdf)) +
    scale_fill_brewer(palette = "Greys") +
    {if(logscale)scale_x_continuous(trans='log10')} +
    ylim(c(0,1)) +
    ylab("") +
    xlab("") +
    annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
    annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
    theme_prj +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")
  p2 <- tidy_success %>%
    filter(timeInt==1) %>%
    mutate(pop_change = 100 * pop_change,
           pop_avg = 100 * pop_avg,
           road_density = road_density/1000) %>%
    rename(value = param) %>%
    distinct(Tag_WMU, Year, value) %>%
    mutate(value = ifelse(rep(logscale,n()),exp(value),value)) %>%
    ggplot(aes(x = value)) +
    geom_density() +
    xlab(xlab) +
    ylab("") +
    {if(logscale)scale_x_continuous(trans='log10')} +
    annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
    annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
    theme_prj +
    theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position="none")
  cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(2,1),
            axis = "rl", align = "v")
}

legend <- cowplot::get_legend(
  tidy_success %>%
    mutate(id = 1:n()) %>%
    ggplot(aes(x = Year, y = pop_change)) +
    ggdist::stat_lineribbon() +
    geom_line(linewidth = 1) +
    scale_fill_brewer(palette = "Greys") +
    labs(fill = "Credible\ninterval", linetype = "")
)

a <- cdf_effects2("pop_change", mod1, mod2, xlab=expression("Change in moose density from average (per 100 km"^2*")"))
b <- cdf_effects2("tagpermooselog_change", mod1, mod2, 
                  xlab="Proportional change in tags per moose from average",
                  logscale = TRUE)
c <- cdf_effects2("precip_change", mod1, mod2, xlab="Change in precipitation from average (mm/month)")
d <- cdf_effects2("pop_avg", mod1, mod2, xlab=expression("Average moose density (per 100 km"^2*")"))
e <- cdf_effects2("tagpermooselog_avg", mod1, mod2, 
                  xlab="Average tags per moose",
                  logscale = TRUE)
g <- cdf_effects2("snow_avg", mod1, mod2, xlab="Average Snow Depth (cm)")
h <- cdf_effects2("road_density", mod1, mod2, xlab="Road Density (km/km^2)")
i <- cdf_effects2("snow_change", mod1, mod2, xlab="Change in Snow Depth from Average (cm)")
j <- cdf_effects2("dhunters_change", mod1, mod2, xlab="Change in hunter density from Average")
k <- cdf_effects2("dhunters_avg", mod1, mod2, xlab="Average hunter density")


p1 <- ggpubr::ggarrange(
  d + ylab("") + ggtitle("Between-WMU effects") + theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill = "Credible interval"),
  a + ylab("") + ggtitle("Within-WMU effects") + theme(plot.title = element_text(hjust = 0.5)),
  e + ylab(""),
  b + ylab(""),
  g + ylab(""),
  i + ylab(""),
  labels = c("(A)","","(B)","","(C)",""),
  vjust = 1.2,
  nrow = 3, ncol = 2, 
  legend = "right",
  align = "hv")

ggpubr::annotate_figure(ggpubr::ggarrange(p1,legend,widths = c(2,0.3)),
                left = "Harvest success rate")
ggsave("cond_effects_main.jpg", height = 8, width = 9, bg = "white")


# Continuous Season Length ----
temp <- tidy_success %>%
  group_by(Tag_Moose, Tag_Res, Tag_FA) %>%
  select(-timeInt, -timeInt_poly, -timeInt_poly1, -timeInt_poly2, -n, -y) %>%
  distinct() %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = round(modelr::seq_range(c(5,90), 20)),
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    timeInt = 1:max(Season_Length),
    Year = "1",
    n = 1, 
    Tag_WMU = "13"
  ) %>%
  filter(Tag_Moose == "Bull",
         Tag_FA == "G",
         Tag_Res == "Res",
         timeInt <= Season_Length) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
    timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
    timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3],
    pdays = timeInt/Season_Length
  ) %>%
  tidybayes::add_epred_rvars(mod1, allow_new_levels = TRUE, re_formula = NA) %>%
  rename(.epred1 = .epred) %>%
  tidybayes::add_epred_rvars(mod2, allow_new_levels = TRUE, re_formula = NA) %>%
  rename(.epred2 = .epred) %>%
  group_by(Tag_WMU, Tag_Moose, Tag_Res, Tag_FA, Season_Length) %>%
  arrange(timeInt) %>%
  mutate(n = 1000) %>%
  summarise(
    cdf = rvar_func(timeInt, .epred1, .epred2),
  ) %>%
  group_by(Tag_Moose, Tag_Res, Tag_FA, Season_Length) %>%
  summarise(
    cdf = rvar_mean(cdf)
  ) %>%
  filter(Season_Length <= 60 | Tag_Res == "Res") %>%
  ggplot(aes(x = Season_Length)) +
  ggdist::stat_dist_lineribbon(aes(dist = cdf)) +
  scale_fill_brewer(palette = "Greys") +
  ylim(c(0,1)) +
  ylab("Harvest success rate") +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

p2 <- tidy_success %>%
  filter(timeInt==1) %>%
  group_by(Tag_WMU, Tag_FA, Tag_Res, Season_Length) %>%
  mutate(tag = ifelse(Tag_Res == "Res","Resident tag", "Outfitter tag"),
         tag = factor(tag, levels = c("Resident tag", "Outfitter tag"))) %>%
  ggplot(aes(x = Season_Length)) +
  geom_histogram() +
  scale_fill_manual(values = c("black", "darkgrey")) +
  xlab("Season length (days)") +
  ylab("") +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
cowplot::plot_grid(temp +
                     labs(fill = "Credible\ninterval"),
                   p2, ncol = 1, rel_heights = c(3,1),
                   axis = "rl", align = "v")

ggsave("season_length.jpg",width = 6, height = 4)


# Between tag types ----
tidy_success %>%
  dplyr::select(-timeInt, -timeInt_poly, -timeInt_poly1, -timeInt_poly2, -n, -y) %>%
  distinct() %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = 40,
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    Year = "1",
    n = 1, 
    Tag_WMU = "13",
    Tag_Moose = c("Bull", "Cow"),
    Tag_Res,
    Tag_FA,
    timeInt = 1:Season_Length
  ) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
    timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
    timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3],
    pdays = timeInt/Season_Length
  ) %>%
  tidybayes::add_epred_rvars(mod1, allow_new_levels = TRUE, re_formula = NA) %>%
  rename(.epred1 = .epred) %>%
  tidybayes::add_epred_rvars(mod2, allow_new_levels = TRUE, re_formula = NA) %>%
  rename(.epred2 = .epred) %>%
  group_by(Tag_WMU, Tag_Moose, Tag_Res, Tag_FA) %>%
  arrange(timeInt) %>%
  summarise(
    cdf = rvar_func(timeInt, .epred1, .epred2),
  ) %>%
  mutate(
    Moose = str_c(Tag_Moose, " tag"),
    Tag_Res = fct_recode(Tag_Res, "Resident\nlottery tag" = "Res", "Tourist\noutfitter tag" = "Tour"),
    Firearm = fct_recode(Tag_FA, "Gun tag" = "G", "Bow tag" = "A")
  ) %>%
  ggplot(aes(x = Tag_Res, fill = Firearm)) +
  ggdist::stat_dist_pointinterval(aes(dist = cdf), position = "dodge", 
                          colour = "black", shape = 21, 
                          size = 6, stroke = 0.5) +
  ylim(c(0,1)) +
  ylab("Harvest success rate") +
  xlab("") +
  facet_wrap(~Moose) +
  labs(fill = "") +
  scale_fill_manual(values = c("white", "black")) +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj
ggsave("condeffects_tagtype.jpg",width = 5,height=4)



# Hazard and Survival Functions ----
s1 <- tidy_success %>%
  group_by(Tag_Moose, Tag_Res, Tag_FA) %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = c(5,30,60),
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    Year = "1",
    timeInt = 1:60,
    n = 1, 
    Tag_WMU = "15B" 
  ) %>%
  filter(Tag_Moose == "Bull",
         Tag_FA == "G",
         Tag_Res == "Res",
         timeInt <= Season_Length) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
    timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
    timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3]
  ) %>%
  tidybayes::add_epred_rvars(mod1, ndraws = 3000, allow_new_levels = TRUE, 
                  re_formula = NA) %>%
  group_by(Tag_WMU, Tag_Res, Tag_Moose, Tag_FA, Year, Season_Length) %>%
  arrange(timeInt) %>%
  summarise(
    cdf = rvar_mult(timeInt, .epred/n),
    timeInt = as.list(timeInt),
  ) %>%
  unnest(c(cdf, timeInt)) %>%
  group_by(Tag_Res, Tag_Moose, Tag_FA, timeInt, Season_Length) %>%
  summarise(cdf = rvar_mean(cdf)) %>%
  mutate(Season_Length = ifelse(Season_Length == 5,
                                "5-day season",
                                ifelse(Season_Length == 30,
                                       "30-day season",
                                       "60-day season")),
         Season_Length = factor(Season_Length,
                                levels = c("5-day season",
                                           "30-day season",
                                           "60-day season"))) %>%  
  ggplot(aes(x = timeInt)) +
  ggdist::stat_dist_lineribbon(aes(dist = cdf)) +
  geom_point(aes(alpha = timeInt %in% seq(0,60,5), y = median(cdf)), 
             colour = "black", size = 2) +
  scale_fill_brewer(palette = "Greys") +
  scale_alpha_manual(values = c(0, 1)) +
  ylab("Cumulative success probability") +
  xlab("Days") +
  labs(fill = "Credible\ninterval")+
  facet_wrap(~Season_Length, scales = "free_x") +
  ylim(c(0,1)) +
  guides(alpha = "none") +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj

h1 <- tidy_success %>%
  group_by(Tag_Moose, Tag_Res, Tag_FA) %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = c(5,30,60),
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    Year = "1",
    timeInt = 1:60,
    n = 1, 
    Tag_WMU = "15B" 
  ) %>%
  filter(Tag_Moose == "Bull",
         Tag_FA == "G",
         Tag_Res == "Res",
         timeInt <= Season_Length) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
    timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
    timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3]
  ) %>%
  tidybayes::add_epred_rvars(mod1, ndraws = 3000, allow_new_levels = TRUE, 
                  re_formula = NA) %>%
  group_by(Tag_Res, Tag_Moose, Tag_FA, timeInt, Season_Length) %>%
  summarise(cdf = rvar_mean(.epred)) %>%
  mutate(Season_Length = ifelse(Season_Length == 5,
                                "5-day season",
                                ifelse(Season_Length == 30,
                                       "30-day season",
                                       "60-day season")),
         Season_Length = factor(Season_Length,
                                levels = c("5-day season",
                                           "30-day season",
                                           "60-day season"))) %>%  
  ggplot(aes(x = timeInt)) +
  ggdist::stat_dist_lineribbon(aes(dist = cdf)) +
  geom_point(aes(alpha = timeInt %in% seq(0,60,5), y = median(cdf)), 
             colour = "black", size = 2) +
  scale_fill_brewer(palette = "Greys") +
  scale_alpha_manual(values = c(0, 1)) +
  ylab("Conditional success probability") +
  xlab("Day of hunt") +
  labs(fill = "Credible\ninterval")+
  facet_wrap(~Season_Length, scales = "free_x") +
  guides(alpha = "none") +
  ylim(c(0,1)) +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj

s2 <- tidy_unsuccess %>%
  group_by(Tag_Moose, Tag_Res, Tag_FA) %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = c(5,30,60),
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    Year = "1",
    timeInt = 1:60,
    n = 1, 
    Tag_WMU = "15B" 
  ) %>%
  filter(Tag_Moose == "Bull",
         Tag_FA == "G",
         Tag_Res == "Res",
         timeInt < Season_Length) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
    timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
    timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3],
    pdays = timeInt/Season_Length
  ) %>%
  tidybayes::add_epred_rvars(mod2, ndraws = 3000, allow_new_levels = TRUE, 
                  re_formula = NA) %>%
  group_by(Tag_WMU, Tag_Res, Tag_Moose, Tag_FA, Year, Season_Length) %>%
  arrange(timeInt) %>%
  summarise(
    cdf = rvar_mult(timeInt, .epred/n),
    timeInt = as.list(timeInt),
  ) %>%
  unnest(c(cdf, timeInt)) %>%
  group_by(Tag_Res, Tag_Moose, Tag_FA, timeInt, Season_Length) %>%
  summarise(cdf = rvar_mean(cdf)) %>%
  mutate(Season_Length = ifelse(Season_Length == 5,
                                "5-day season",
                                ifelse(Season_Length == 30,
                                       "30-day season",
                                       "60-day season")),
         Season_Length = factor(Season_Length,
                                levels = c("5-day season",
                                           "30-day season",
                                           "60-day season"))) %>%  
  ggplot(aes(x = timeInt)) +
  ggdist::stat_dist_lineribbon(aes(dist = cdf)) +
  geom_point(aes(alpha = timeInt %in% seq(0,60,5), y = median(cdf)), 
             colour = "black", size = 2) +
  scale_fill_brewer(palette = "Greys") +
  scale_alpha_manual(values = c(0, 1)) +
  ylab("Cumulative probability of ending hunt") +
  xlab("Days") +
  labs(fill = "Credible\ninterval")+
  facet_wrap(~Season_Length, scales = "free_x") +
  ylim(c(0,1)) +
  guides(alpha = "none") +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj

h2 <- tidy_unsuccess %>%
  group_by(Tag_Moose, Tag_Res, Tag_FA) %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = c(5,30,60),
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    Year = "1",
    timeInt = 1:60,
    n = 1, 
    Tag_WMU = "15B" 
  ) %>%
  filter(Tag_Moose == "Bull",
         Tag_FA == "G",
         Tag_Res == "Res",
         timeInt < Season_Length) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_log_poly1 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,1],
    timeInt_log_poly2 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,2],
    timeInt_log_poly3 = predict(tidy_success$timeInt_log_poly, log(timeInt))[,3],
    pdays = timeInt/Season_Length
  ) %>%
  tidybayes::add_epred_rvars(mod2, ndraws = 3000, allow_new_levels = TRUE, 
                  re_formula = NA) %>%
  group_by(Tag_Res, Tag_Moose, Tag_FA, timeInt, Season_Length) %>%
  summarise(cdf = rvar_mean(.epred)) %>%
  mutate(Season_Length = ifelse(Season_Length == 5,
                                "5-day season",
                                ifelse(Season_Length == 30,
                                      "30-day season",
                                      "60-day season")),
         Season_Length = factor(Season_Length,
                                levels = c("5-day season",
                                           "30-day season",
                                           "60-day season"))) %>%  
  ggplot(aes(x = timeInt)) +
  ggdist::stat_dist_lineribbon(aes(dist = cdf)) +
  geom_point(aes(alpha = timeInt %in% seq(0,60,5), y = median(cdf)), 
             colour = "black", size = 2) +
  scale_fill_brewer(palette = "Greys") +
  scale_alpha_manual(values = c(0, 1)) +
  ylab("Conditional probability of ending hunt") +
  xlab("Day of hunt") +
  labs(fill = "Credible\ninterval")+
  facet_wrap(~Season_Length, scales = "free_x") +
  guides(alpha = "none") +
  ylim(c(0,1)) +
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  annotation_custom(grid::linesGrob(x = c(0, 0), gp = grid::gpar(lwd = 3))) +
  theme_prj


ggpubr::ggarrange(
  h1,
  h2,
  s1,
  s2,
  nrow = 2, ncol = 2, 
  labels = c("(A)","(B)","(C)","(D)"),
  common.legend = TRUE,
  legend = "right",
  align = "hv")
ggsave("hazardsurvivalfunction.jpg", width = 8, height = 7)

