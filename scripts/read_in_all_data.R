###########################
###########################
## This script will pull the data into the working environment for those with 
## access to the data repo for the DFO St. Lawrence data
###########################
###########################
## date: 2020-11-23
## author: Cole Brookson
###########################
###########################

require(here)
require(stringr)
require(tidyverse)
library(patchwork)

#set the directory using this super hacky method to keep the data private
dir = paste0(str_sub(here(), 1, -9),'ciee-stlawrence/data')
setwd(dir)

#load all data into the environment
explan = readRDS("explan.RDS")
number = readRDS("number.RDS")
gulf = readRDS("gulf.RDS")
catch_total = readRDS("catchTotal.RDS")
fish_rich = readRDS("fishRich.RDS")
dem = readRDS("dem.RDS")
invert_rich = readRDS("invertRich.RDS")
strata = readRDS("strata.RDS")
weight = readRDS("weight.RDS")
species_names = read_csv('species_names.csv')

rm(dir)

#round the numbers of the abundance
# number@data = round(number@data)
# for(i in 1:ncol(number@data)) {
#   if(class(number@data[,i]) == "numeric"){
#     number@data[,i] = as.integer(number@data[,i])
#   }
# }

