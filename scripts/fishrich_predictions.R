###########################
###########################
## This script will try to generate predictions from simple models using  
## INLA for fish richness for the DFO St. Lawrence dataset
###########################
###########################
## date: 2020-12-08
## author: Cole B. Brookson
###########################
###########################

# The hypotheses models are:

# 1) btemp (energy hypothesis)
# 2) ice cover, bSalinity (productivity hypothesis)
# 3) sTemp, bTemp (climate stability)
# 4) slope, bSalinity (spatial heterogeneity hypothesis)
# 5) depth (stress hypothesis)

# full model: fishrich ~ bTemp + sTemp + ice cover + bSalinity + slope + depth

# set-up =======================================================================

library(INLA)
library(tidyverse)
library(ggregplot) 
library(inlabru)

# read in model outputs from model runs done in `fishrich.R`
m_full = readRDS(here("./model-output/fishrich_full.rds"))
m_energy = readRDS(here("./model-output/fishrich_energy.rds"))
m_prod = readRDS(here("./model-output/fishrich_prod.rds"))
m_clim = readRDS(here("./model-output/fishrich_clim.rds"))
m_sphet = readRDS(here("./model-output/fishrich_sphet.rds"))
m_stress = readRDS(here("./model-output/fishrich_stress.rds"))







