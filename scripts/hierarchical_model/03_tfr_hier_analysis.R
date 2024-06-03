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

### load WMU shapefile
WMU_shp <- readRDS("data/processed/wmu_shp_simp.rds")

## these calculate model r squared (no way of taking out the effect of group 
### level effects to see the explanatory power of covariates alone)
performance::r2_bayes(mod1)
performance::r2_bayes(mod2)
brms::bayes_R2(mod1, re_formula = NA)
brms::bayes_R2(mod2, re_formula = NA)
brms::pp_check(mod1, type = "bars", ndraws = 1000) + xlim(-2,10)
brms::pp_check(mod2, type = "bars", ndraws = 1000) + xlim(-2,10)

temp <- tidy_unsuccess %>%
  filter(timeInt <= 20)

temp %>%
  mutate(Tag_Moose = droplevels(Tag_Moose),
         Tag_WMU = droplevels(Tag_WMU),
         Year = droplevels(Year),
         pdays = timeInt/Season_Length) %>%
  add_predicted_draws(mod2, ndraws=100) %>%
  group_by(.draw, .prediction, timeInt) %>%
  summarise(n_rep = n()) %>%
  ggplot(aes(x = .prediction)) +
  geom_histogram(data = temp, aes(x = y), binwidth = 1) +
  stat_pointinterval(aes(y = n_rep),.width = c(.50, .80, .95, .99)) +
  facet_wrap(~timeInt, scales = "free_y") +
  xlim(-2,10) 

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


loo1 <- loo::loo(mod1)
loo2 <- loo::loo(mod2)
### p=49 or p=490 or somewhere in between
loo_r21 <- loo_R2(mod1)
loo_r22 <- loo_R2(mod2)
saveRDS(loo_r21, "output/hierarchical_model_looR2_mod1.rds")
saveRDS(loo_r22, "output/hierarchical_model_looR2_mod2.rds")
loo_r21 <- readRDS("output/hierarchical_model_looR2_mod1.rds")
loo_r22 <- readRDS("output/hierarchical_model_looR2_mod2.rds")
## this does a different method
performance::r2_nakagawa(mod1, moment)
rstanarm_bayes_r2_fix <- function(object, re.form = NULL) {
  mu_pred <- rstanarm::posterior_predict(object, re.form = re.form)
  y <- insight::get_response(object)[,1]
  rstantools::bayes_R2(mu_pred, y = y)
}
median(rstanarm_bayes_r2_fix(mod1, re.form = NA))
median(rstanarm_bayes_r2_fix(mod2, re.form = NA))

# Build coefficient table ----
### specify parameter names
coef_levels <- c("Season length",
                 "Road density",
                 "Moose density (temporal)","Moose density (spatial)", 
                 "Tags per moose (temporal)","Tags per moose (spatial)",
                 "Hunter density (temporal)", "Hunter density (spatial)",
                 "Snow depth (temporal)", "Snow depth (spatial)",
                 "Precipitation (temporal)","Temperature (temporal)")
coef_levels2 <- c("Intercept",
                  "Tag sex (cow)",
                  "Tag residency (tourism)",
                  "Tag firearm (bow)",
                  "Tag sex:residency",
                  "Tag sex:firearm",
                  "Tag residency:firearm",
                  "Tag sex:residency:firearm",
                  "Day of hunt (linear term)",
                  "Day of hunt (quadratic term)",
                  "Day of hunt/Season length")
coef_levels3 <- c("sd CAR",
                  "sd WMU intercept",
                  "sd WMU tag sex (cow)",
                  "sd WMU tag residency (tourism)",
                  "sd WMU tag firearm (bow)",
                  "sd WMU tag sex:residency",
                  "sd WMU tag sex:firearm",
                  "sd WMU tag residency:firearm",
                  "sd WMU tag sex:residency:firearm",
                  "sd Year intercept")
coef_levels4 <- c("Day of hunt (linear term)",
                  "Day of hunt (quadratic term)",
                  "Day of hunt/Season length",
                  "Resident Gun Bull",
                  "Resident Gun Cow",
                  "Resident Bow Bull",
                  "Tourism Gun Bull",
                  "Resident Bow Cow",
                  "Tourism Bow Bull",
                  "Tourism Gun Cow",
                  "Tourism Bow Cow")
colour <- c("Spatial","Spatial","Temporal","Spatial",
            "Temporal","Spatial","Temporal","Spatial",
            "Temporal","Spatial","Temporal")


### example coefficient plot to take ggplot theme from
coef_plot <- brms::mcmc_plot(mod1, 
                       type = "intervals",
                       variable = c('b_scaleSeason_Length', 
                                    'b_scaleroad_density',
                                    'b_scalepop_change', 'b_scalepop_avg', 
                                    'b_scaletagpermooselog_change', 'b_scaletagpermooselog_avg', 
                                    'b_scaledhunters_change', 'b_scaledhunters_avg', 
                                    'b_scalesnow_change', 'b_scalesnow_avg',
                                    'b_scaleprecip_change','b_scaletemp_change')
) 

### summarise parameter estimates for beta coefficients for both models
beta1 <- mod1 %>%
  tidybayes::spread_draws(`b_scaleSeason_Length`, 
               `b_scaleroad_density`,
               `b_scalepop_change`,  `b_scalepop_avg`, 
               `b_scaletagpermooselog_change`, `b_scaletagpermooselog_avg`, 
               `b_scaledhunters_change`, `b_scaledhunters_avg`, 
               `b_scalesnow_change`, `b_scalesnow_avg`,
               `b_scaleprecip_change`,`b_scaletemp_change`
  ) %>%
  tidybayes::summarise_draws(~quantile(.x, probs = c(0.025,0.25,0.5,0.75,0.975)),
                  default_convergence_measures()) %>%
  rename(parameter = variable,
         m = `50%`,
         l = `25%`,
         h = `75%`,
         ll = `2.5%`,
         hh = `97.5%`) %>%
  mutate(parameter = factor(
    c("Season length",
      "Road density",
      "Moose density (temporal)","Moose density (spatial)", 
      "Tags per moose (temporal)","Tags per moose (spatial)",
      "Hunter density (temporal)", "Hunter density (spatial)",
      "Snow depth (temporal)", "Snow depth (spatial)",
      "Precipitation (temporal)","Temperature (temporal)"), 
    levels = rev(coef_levels)
  ),
  model = "Probability of success"
  )
beta2 <- mod2 %>%
  tidybayes::spread_draws(`b_scaleSeason_Length`, 
               `b_scaleroad_density`,
               `b_scalepop_change`, `b_scalepop_avg`, 
               `b_scaledhunters_change`, `b_scaledhunters_avg`, 
               `b_scalesnow_change`,`b_scalesnow_avg`,
               `b_scaleprecip_change`,`b_scaletemp_change`
  ) %>%
  tidybayes::summarise_draws(~quantile(.x, probs = c(0.025,0.25,0.5,0.75,0.975)),
                  default_convergence_measures()) %>%
  rename(parameter = variable,
         m = `50%`,
         l = `25%`,
         h = `75%`,
         ll = `2.5%`,
         hh = `97.5%`) %>%
  mutate(parameter = factor(
    c("Season length",
      "Road density",
      "Moose density (temporal)", "Moose density (spatial)", 
      "Hunter density (temporal)", "Hunter density (spatial)",
      "Snow depth (temporal)", "Snow depth (spatial)",
      "Precipitation (temporal)","Temperature (temporal)"), 
    levels = rev(coef_levels)
  ),
  model = "Probability of ending hunt"
  )

### summarise parameter estimates for categorical and day of hunt effects 
### for both models
inter1 <- mod1 %>%
  tidybayes::spread_draws(b_Intercept, b_Tag_MooseCow,b_Tag_ResTour,
               b_Tag_FAA,`b_Tag_ResTour:Tag_MooseCow`,
               `b_Tag_MooseCow:Tag_FAA`,
               `b_Tag_ResTour:Tag_FAA`,
               `b_Tag_ResTour:Tag_MooseCow:Tag_FAA`,
               b_scaletimeInt_poly1,
               b_scaletimeInt_poly2
  ) %>%
  tidybayes::summarise_draws(~quantile(.x, probs = c(0.025,0.25,0.5,0.75,0.975)),
                  default_convergence_measures()) %>%
  rename(parameter = variable,
         m = `50%`,
         l = `25%`,
         h = `75%`,
         ll = `2.5%`,
         hh = `97.5%`) %>%
  mutate(parameter = factor(
    c("Intercept",
      "Tag sex (cow)",
      "Tag residency (tourism)",
      "Tag firearm (bow)",
      "Tag sex:residency",
      "Tag sex:firearm",
      "Tag residency:firearm",
      "Tag sex:residency:firearm",
      "Day of hunt (linear term)",
      "Day of hunt (quadratic term)"), 
    levels = rev(coef_levels2)
  ),
  model = "Probability of success"
  )

inter2 <- mod2 %>%
  tidybayes::spread_draws(b_Intercept, b_Tag_MooseCow,b_Tag_ResTour,
               b_Tag_FAA,`b_Tag_ResTour:Tag_MooseCow`,
               `b_Tag_MooseCow:Tag_FAA`,
               `b_Tag_ResTour:Tag_FAA`,
               `b_Tag_ResTour:Tag_MooseCow:Tag_FAA`,
               b_scaletimeInt_poly1,
               b_scaletimeInt_poly2,
               b_scalepdays
  )  %>%
  tidybayes::summarise_draws(~quantile(.x, probs = c(0.025,0.25,0.5,0.75,0.975)),
                  default_convergence_measures()) %>%
  rename(parameter = variable,
         m = `50%`,
         l = `25%`,
         h = `75%`,
         ll = `2.5%`,
         hh = `97.5%`) %>%
  mutate(parameter = factor(
    c("Intercept",
      "Tag sex (cow)",
      "Tag residency (tourism)",
      "Tag firearm (bow)",
      "Tag sex:residency",
      "Tag sex:firearm",
      "Tag residency:firearm",
      "Tag sex:residency:firearm",
      "Day of hunt (linear term)",
      "Day of hunt (quadratic term)",
      "Day of hunt/Season length"), 
    levels = rev(coef_levels2)
  ),
  model = "Probability of ending hunt"
  )

### summarise parameter estimates for group-level effects for both models
groups1 <- mod1 %>%
  tidybayes::spread_draws(sdcar,sd_Tag_WMU__Intercept, sd_Tag_WMU__Tag_ResTour,
               sd_Tag_WMU__Tag_MooseCow, sd_Tag_WMU__Tag_FAA,
               `sd_Tag_WMU__Tag_ResTour:Tag_MooseCow`,
               `sd_Tag_WMU__Tag_MooseCow:Tag_FAA`,
               `sd_Tag_WMU__Tag_ResTour:Tag_FAA`,
               `sd_Tag_WMU__Tag_ResTour:Tag_MooseCow:Tag_FAA`,
               `sd_Year__Intercept`
  ) %>%
  tidybayes::summarise_draws(~quantile(.x, probs = c(0.025,0.25,0.5,0.75,0.975)),
                  default_convergence_measures()) %>%
  rename(parameter = variable,
         m = `50%`,
         l = `25%`,
         h = `75%`,
         ll = `2.5%`,
         hh = `97.5%`) %>%
  mutate(parameter = factor(
    c("sd CAR",
      "sd WMU intercept",
      "sd WMU tag sex (cow)",
      "sd WMU tag residency (tourism)",
      "sd WMU tag firearm (bow)",
      "sd WMU tag sex:residency",
      "sd WMU tag sex:firearm",
      "sd WMU tag residency:firearm",
      "sd WMU tag sex:residency:firearm",
      "sd Year intercept"), 
    levels = rev(coef_levels3)
  ),
  model = "Probability of success"
  )

groups2 <- mod2 %>%
  tidybayes::spread_draws(sd_Tag_WMU__Intercept, sd_Tag_WMU__Tag_ResTour,
               sd_Tag_WMU__Tag_MooseCow, sd_Tag_WMU__Tag_FAA,
               `sd_Tag_WMU__Tag_ResTour:Tag_MooseCow`,
               `sd_Tag_WMU__Tag_MooseCow:Tag_FAA`,
               `sd_Tag_WMU__Tag_ResTour:Tag_FAA`,
               `sd_Tag_WMU__Tag_ResTour:Tag_MooseCow:Tag_FAA`,
               `sd_Year__Intercept`
  ) %>%
  tidybayes::summarise_draws(~quantile(.x, probs = c(0.025,0.25,0.5,0.75,0.975)),
                  default_convergence_measures()) %>%
  rename(parameter = variable,
         m = `50%`,
         l = `25%`,
         h = `75%`,
         ll = `2.5%`,
         hh = `97.5%`) %>%
  mutate(parameter = factor(
    c("sd WMU intercept",
      "sd WMU tag sex (cow)",
      "sd WMU tag residency (tourism)",
      "sd WMU tag firearm (bow)",
      "sd WMU tag sex:residency",
      "sd WMU tag sex:firearm",
      "sd WMU tag residency:firearm",
      "sd WMU tag sex:residency:firearm",
      "sd Year intercept"), 
    levels = rev(coef_levels3)
  ),
  model = "Probability of ending hunt"
  )

### combine all parameter estimates and plot coefficient table
beta1 %>%
  bind_rows(beta2) %>%
  bind_rows(inter1) %>%
  bind_rows(inter2) %>%
  bind_rows(groups1) %>%
  bind_rows(groups2) %>%
  mutate(m = format(round(m,2), nsmall = 2),
         ll = format(round(ll,2), nsmall = 2),
         hh = format(round(hh,2), nsmall = 2),
         `95% CI` = str_c("(",ll,", ",hh,")")) %>%
  dplyr::select(-l,-h, -ll, -hh) %>%
  rename(Median = m,
         Parameter = parameter) %>%
  pivot_wider(names_from = model, values_from = c(Median,`95% CI`,rhat,ess_bulk,ess_tail)) %>%
  write_csv("coef_table.csv")


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

# Between WMUs ----
WMU_shp_simp <- readRDS("data/processed/WMU_shp_simp.rds")

tidy_success %>%
  group_by(Tag_WMU, Tag_Moose, Tag_Res, Tag_FA) %>%
  select(-timeInt, -timeInt_poly, -timeInt_poly1, -timeInt_poly2, -n, -y) %>%
  distinct() %>%
  modelr::data_grid(
    pop_avg = mean(pop_avg),
    tagpermooselog_avg = mean(tagpermooselog_avg),
    pop_change = 0,
    tagpermooselog_change = 0,
    Season_Length = mean(Season_Length),
    snow_avg = mean(snow_avg),
    precip_change = 0,
    temp_change = 0,
    snow_change = 0,
    road_density = mean(road_density, na.rm = TRUE),
    dhunters_avg = mean(dhunters_avg),
    dhunters_change = 0,
    timeInt = 1:mean(Season_Length),
    Year = "1",
    n = 1 
  ) %>%
  mutate(
    Tag_Moose = as.character(Tag_Moose),
    timeInt_poly1 = predict(tidy_success$timeInt_poly, timeInt)[,1],
    timeInt_poly2 = predict(tidy_success$timeInt_poly, timeInt)[,2],
    pdays = timeInt/Season_Length
  ) %>%
  tidybayes::add_epred_rvars(mod1, allow_new_levels = TRUE, ndraws = 1000) %>%
  rename(.epred1 = .epred) %>%
  tidybayes::add_epred_rvars(mod2, allow_new_levels = TRUE, ndraws = 1000) %>%
  rename(.epred2 = .epred) %>%
  group_by(Tag_WMU, Tag_Moose, Tag_Res, Tag_FA) %>%
  arrange(timeInt) %>%
  summarise(
    cdf = mean(rvar_func(timeInt, .epred1, .epred2)),
  ) %>%
  mutate(ecdf = mean(cdf)) %>%
  right_join(WMU_shp_simp %>%
               expand_grid(Tag_Moose = c("Bull", "Cow"),
                           Tag_Res = c("Res", "Tour"),
                           Tag_FA = c("G", "A")),
             by = c("Tag_WMU" = "WMU", 
                    "Tag_Moose", "Tag_Res", "Tag_FA")) %>%
  sf::st_as_sf() %>%
  mutate(Tag_type = str_c(Tag_Res, Tag_FA, Tag_Moose, sep = "")) %>%
  ggplot(aes(fill = ecdf)) +
  geom_sf() +
  labs(fill = "Proportion of\nsuccessful\nhunters") +
  scale_fill_viridis_c(option = "plasma",
                       na.value = "white") +
  facet_wrap(~Tag_type, nrow = 2) +
  theme_prj
ggsave(str_c(picf,"pred_WMU_TFR.png"),  height = 8, width = 12)


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



# Recreate Data (Number of successful hunters) ----
data_rep <- tidy_success %>%
  filter(timeInt == 1) %>%
  mutate(Tag_Moose = as.character(Tag_Moose),
         y = NULL,
         n_val = n,
         n = 1,
         Year_val = Year,
         Year = "1",
         timeInt = NULL,
         timeInt_poly = NULL,
         timeInt_poly1 = NULL,
         timeInt_poly2 = NULL) %>%
  distinct() %>%
  expand_grid(timeInt = 1:90) %>%
  filter(timeInt <= Season_Length) %>%
  mutate(timeInt_poly1 = predict(tidy_success$timeInt_poly, timeInt)[,1],
       timeInt_poly2 = predict(tidy_success$timeInt_poly, timeInt)[,2]) %>%
  filter(Year_val == 2019) %>%
  tidybayes::add_epred_rvars(mod1, ndraws = 3000, allow_new_levels = TRUE)

data_rep2 <- tidy_success %>%
  filter(timeInt == 1) %>%
  mutate(Tag_Moose = as.character(Tag_Moose),
         n_val = n,
         n = 1,
         y = NULL,
         Year_val = Year,
         Year = "1",
         timeInt = NULL,
         timeInt_poly = NULL,
         timeInt_poly1 = NULL,
         timeInt_poly2 = NULL) %>%
  distinct() %>%
  expand_grid(timeInt = 1:90) %>%
  filter(timeInt <= Season_Length) %>%
  mutate(timeInt_poly1 = predict(tidy_success$timeInt_poly, timeInt)[,1],
         timeInt_poly2 = predict(tidy_success$timeInt_poly, timeInt)[,2],
         pdays = timeInt/Season_Length) %>%
  filter(Year_val == 2019) %>%
  tidybayes::add_epred_rvars(mod2, ndraws = 3000, allow_new_levels = TRUE)

## Using prob of censoring
temp <- data_rep %>% 
  mutate(Year = Year_val) %>%
  left_join(
    data_rep2 %>%
      dplyr::select(Year_val, Tag_WMU, Tag_Res, Tag_Moose, Tag_FA, .epred, timeInt) %>%
      rename(Year = Year_val, .epred2 = .epred)) %>%
  group_by(Tag_WMU, Tag_Res, Tag_Moose, Tag_FA, Year) %>%
  arrange(timeInt) %>%
  summarise(y_pred = rvar_func2(timeInt, max(n_val), .epred, .epred2))

tidy_success %>%
  group_by(Tag_WMU, Tag_Res, Tag_Moose, Tag_FA, Year) %>%
  summarise(y = sum(y, na.rm = TRUE),
            n = max(n, na.rm = TRUE)) %>%
  right_join(temp) %>%
  group_by(Tag_WMU, Tag_Moose) %>%
  summarise(y_pred = rvar_sum(y_pred),
            y = sum(y, na.rm = TRUE),
            n2 = sum(n, na.rm = TRUE),
            p = mean(y/n, na.rm = TRUE)) %>%
  ggplot(aes(y = y/n2, x = Tag_WMU, fill = "observed")) +
  ggdist::stat_dist_pointinterval(aes(dist = y_pred/n2, fill = "modelled")) +
  geom_point(shape = 21, size = 2) +
  facet_grid(~Tag_Moose) +
  ylab("Harvest success rate") +
  xlab("WMU") +
  labs(fill = NULL) +
  scale_fill_manual(values = c("black", "red")) +
  theme_prj

ggsave(str_c(picf,"modelvsobs_2019.png"), width = 9, height = 5)
