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
library(raster)
library(fields)

# read in model outputs from model runs done in `fishrich.R`
m_full = readRDS(here("./model-output/fishrich_full.rds"))
m_energy = readRDS(here("./model-output/fishrich_energy.rds"))
m_prod = readRDS(here("./model-output/fishrich_prod.rds"))
m_clim = readRDS(here("./model-output/fishrich_clim.rds"))
m_sphet = readRDS(here("./model-output/fishrich_sphet.rds"))
m_stress = readRDS(here("./model-output/fishrich_stress.rds"))

### NOTE ## -------------------------------------------------------------------
# the productivity hypothesis (m_prod) is the best according to model selection
# so we'll go with that one for now, and try to add a spatial component
## ----------------------------------------------------------------------------

# create the mesh ==============================================================

# make the spatial component
xy = coordinates(explan_ice)
colnames(xy) = c("locx", "locy")

# temporal component
time = as.numeric(explan_ice@endTime)/60 # Number of minutes since Jan 1, 1970 
timeBasis = time - 884844 # Time starting at 1

# build temporal mesh
timeGr = 11
timeSeq = seq(min(timeBasis),
              max(timeBasis),
              length = timeGr)
meshTime = inla.mesh.1d(timeSeq)

#build spatial mesh
maxEdge = 35000
meshSpace = inla.mesh.2d(boundary = gulf, 
                         max.edge=maxEdge * c(0.5,2),
                         cutoff=maxEdge/4,
                         offset = c(10000,20000))

# SPDE and projection matrix ===================================================

# make the stochastic partial differential equation object
SPDE = inla.spde2.pcmatern(mesh=meshSpace,
                           alpha=2,
                           prior.range=c(1000, 0.5),
                           prior.sigma=c(35, 0.5))

# make the priors
hSpec = list(theta=list(prior='pccor1', param=c(0.2, 0.9)))
precPrior = list(prior='pc.prec', param=c(1, 0.01)) 

# create the spatial field 
Field = inla.spde.make.index("field", n.spde=SPDE$n.spde,
                             n.group = meshTime$n)

# create the projection matirx
A = inla.spde.make.A(meshSpace,
                     loc=coordinates(explan_ice@sp),
                     group = timeBasis,
                     group.mesh = meshTime)
Alist = as.list(rep(1,3))
Alist[[1]] = A

# make the effects from the data
effect_data = list(Field,
              bSal = explan_ice@data$bottom.salinity,
              ice = explan_ice@data$ice.duration)

# put all items in a data stack
stack_data = inla.stack(data=list(fish = fish_rich@data$species.fish.number),
                        A = Alist,
                        effects = effect_data,
                        tag = "basis")

# set the predictions ==========================================================

# new data for the predictions
new_data = expand.grid(bSal = seq(min(explan_ice@data$bottom.salinity, 
                                      na.rm = TRUE), 
                                  max(explan_ice@data$bottom.salinity,
                                      na.rm = TRUE), 
                                  length.out = 10),
                       ice = seq(min(explan_ice@data$ice.duration),
                                 max(explan_ice@data$ice.duration),
                                 length.out = 10))

# make effect for predictions
effect_pred = list(bSal = explan_ice@data$bottom.salinity,
                   ice = explan_ice@data$ice.duration)


# make stack for predictions
stack_pred = inla.stack(data=list(fish = NA),
                        A = as.list(rep(1,2)),
                        effects = effect_pred,
                        tag = "pred_fix")

# make an empty raster for the prediction coordinates
raster_empty = raster(xmn = min(xy[,1]), xmx = max(xy[,1]),
                      ymn = min(xy[,2]), ymx = max(xy[,2]),
                      resolution = 25)







# productivity hypothesis ----
prod <- fish ~ 0 + bSal + ice
m_prod <- quick_inla(prod)

