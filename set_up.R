######################################################
############ Moose Tag Fill Rates Project ############ 
######################################################
### Understanding what factors impact hunter success rates and developing
### a model to predict success rates for upcoming hunting seasons

######################################################
### Script to configure settings and install required packages

### install required packages
pkgs <- function(){
  pkgs <- c("tidyverse","sf","discSurv","lubridate","spdep",
            "zoo","plotly","gstat","loo","brms","tidybayes",
            "posterior")
  not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  install.packages(not_installed)
}
pkgs()

library(tidyverse)

### ggplot theme
theme_prj <-theme_bw(10) + theme(panel.grid = element_blank(),
                                 panel.border = element_blank(),
                                 strip.background = element_blank())

### factor levels ----
moose_level <- c("Bull", "Cow", "Calf")
firearm_level <- c("A", "G")
residence_level <- c("Res", "Tour")
logical_level <- c("1", "2")
logical_level2 <- c("Y", "N")
sex_level <- c("Male", "Female")
age_level <- c("Adult", "Calf")

WMU_level <- c('1A','1B', '1C', '1D', '2', '3', '4', '5', '6', '7A', '7B', '8',
                     '9A', '9B', '10', '11A', '11B', '11C', '12A', '12B', '13', '14', 
                     '15A', '15B', '16A', '16B', '16C', '17', '18A', '18B', '19',  
                     '20', '21A', '21B', '22', '23', '24', '25', 
                     '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', 
                     '36', '37', '38', '39', '40', '41', '42', '43A', '43B', '44', 
                     '45', '46', '47', '48', '49', '50', '51', '53', '54', 
                     '55A', '55B', '56', '57', '58', '59', '60', '61', '62', '63', 
                     '64A', '64B', '65', '66A', '66B', '67', '68A', '68B', '69A-1', 
                     '69A-2', '69A-3', '69B', '70', '71', '72A', '72B', '73', '74A', '74B', 
                     '75', '76A', '76B', '76C', '76D', '76E', '77A', '77B', '77C', '78A', 
                     '78B', '78C', '78D', '78E', '79C', '79D', '80', '81A', 
                     '81B', '82A', '82B', '82C', '83A', '83B', '83C', '84', '85A', '85B', 
                     '85C', '86A', '86B', '87A', '87B', '87C', '87D', '87E', '88', '89A', 
                     '89B', '90A', '90B', '91A', '91B', '92A', '92B', '92C', '92D', '93A', 
                     '93B', '93C', '94A', '94B', '95')
wday_level <- c("Sun", "Mon", "Tue", "Wed", 
                "Thu", "Fri", "Sat")
wday_level_cor <- c("1" = "Sun", "2" = "Mon", "3" = "Tue", "4" = "Wed", 
                    "5" = "Thu", "6" = "Fri", "7" = "Sat")
month_level <- c("Sep", "Oct", "Nov", "Dec")
