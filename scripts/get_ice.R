###########################
###########################
## This script will pull the duration of the ice cover data from 
## https://github.com/duplisea/gslea and manipulate it to match with the explan
## dataset for the DFO St. Lawrence data
###########################
###########################
## date: 2020-12-08
## author: Katherine Hébert
###########################
###########################

# load packages
# devtools::install_github("duplisea/gslea", build_vignettes = TRUE)
library(gslea)
library(tidyverse)
library(here)

# load explanatory variables dataset
source(here('./scripts/read_in_all_data.R'))

# Duration of the ice season (number of days)
EA.plot.f(variables = find.vars.f("ice.duration"), years = 1800:2100, EARs=5:6, pch = 16)
ice <- EA.query.f(find.vars.f("ice.duration"), years = 1800:2100, EARs=5:6)
ice$variable <- as.character(ice$variable)
ice <- rename(ice, "ice.duration" = "value")

# add year column to explan for matching by year ----

# reassigning this to make the code a bit simpler (I am lazy!)
df <- explan@data
df <- mutate(df, year = lubridate::year(explan@endTime))

# add EAR column to match with ice cover ----

# note: EAR = Regions 5 and 6 covers the entire sGSL
# (matched EAR to strata by comparing the "Study area without names.tiff" and the 
# map on the gslea Github. Not sure if the match-ups are correct/ if this is the
# best way of doing this...)
df$EAR <- 5
# only these 3 strata seem to be in the 6 EAR:
df$EAR[which(df$stratum %in% c(421, 402, 432))] <- 6

# match explan and glsea dataset
explan@data <- left_join(df, ice, by = c("year", "EAR")) %>% subset(select = -c(variable))

# save
#set the directory using this super hacky method to keep the data private
dir = paste0(str_sub(here(), 1, -9),'ciee-stlawrence/data')
setwd(dir)
saveRDS(explan, "explan_ice.RDS")
