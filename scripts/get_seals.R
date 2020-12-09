###########################
###########################
## This script will add estimated number of seal-years (thousands) spent in the 
# southern Gulf of St. Lawrence in https://cdnsciencepub.com/doi/pdf/10.1139/cjfas-2017-0190
## and manipulate it to match with the explan dataset for the DFO St. Lawrence data
###########################
###########################
## date: 2020-12-09
## author: Katherine Hébert
###########################
###########################

# load packages
library(tidyverse)

# load datasets
explan <- readRDS("data/explan_ice.RDS") # this already has a year column
seals <- read_csv("data/seal-sGSL.csv")

# rename variable
seals <- rename(seals, "sealyears" = "Gulf")

# match explan and seals dataset
explan@data <- left_join(explan@data, seals, by = c("year"))

# save
saveRDS(explan, "data/explan_ice_seals.RDS")
