###########################
###########################
## This script will work on the model fitting process for redfish
## in the DFO St. Lawrence data
###########################
###########################
## date: 2020-11-13
## author: Cole Brookson
###########################
###########################

#libraries
library(INLA)
library(here)
library(animation)
library(raster)

#read in all data
source(here('./scripts/read_in_all_data.R'))

#make some meshes
redfish = number@data$`REDFISH UNSEPARATED`

maxEdge = 25000
meshSpace = inla.mesh.2d(boundary = gulf, 
                         max.edge = maxEdge*c(0.5,2),
                         offset = c(100000,20000),
                         cutoff = maxEdge/2)
par(mar = c(1,1,1,1))
plot(meshSpace, main = "", asp = 1)

meshSpace$n

time = as.numeric(fish_rich@endTime)/60
timeBasis = time - 884844

timeGr = 11
timeSeq = seq(min(timeBasis), max(timeBasis), length = timeGr)

meshTime = inla.mesh.1d(timeSeq)
meshSpace$n * meshTime$n


SPDE = inla.spde2.pcmatern(mesh=meshSpace,
                           alpha=2,
                           prior.range=c(100, 0.5),
                           prior.sigma=c(10, 0.5))


hSpec = list(theta=list(prior='pccor1', param=c(0.2, 0.9)))
precPrior = list(prior='pc.prec', param=c(1, 0.01)) 

Field = inla.spde.make.index("field", n.spde=SPDE$n.spde,
                             n.group = meshTime$n)


# For estimation
A = inla.spde.make.A(meshSpace,
                     loc=coordinates(number@sp),
                     group = timeBasis,
                     group.mesh = meshTime)

Alist = as.list(rep(1,3))
Alist[[1]] = A

effect = list(Field,
              bTemp =  explan@data$bottom.temperature,
              bSal = explan@data$bottom.salinity)

Stack_redfish = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                   A = Alist,
                   effects = effect,
                   tag="basis")


form_redfish = redfish ~ 0 + bTemp + bSal +
  f(field, model=SPDE,
    group = field.group,
    control.group=list(model='ar1', hyper=hSpec))

model_redfish = inla(form_redfish,
             data = inla.stack.data(Stack_redfish),
             family="gaussian",
             control.family =list(link="identity", hyper = list(theta=precPrior)),
             control.predictor=list(A=inla.stack.A(Stack_redfish),
                                    compute=TRUE, link = 1),
             control.compute=list(waic=TRUE),
             verbose = FALSE)

summary(model_redfish)

# Dimention of the raster
stepsize = 1000
rangeX = range(meshSpace$loc[,1])
rangeY = range(meshSpace$loc[,2])
nxy = round(c(diff(rangeX), 
               diff(rangeY)) / stepsize)

# Define basis of the map
mapBasis = inla.mesh.projector(meshSpace,
                                xlim = rangeX,
                                ylim = rangeY,
                                crs = crs(gulf))

# Calculate prediction
mapMean = vector("list", length = timeGr)
map.025 = vector("list", length = timeGr)
map.975 = vector("list", length = timeGr)

for(i in 1:timeGr){
  # Model prediction
  fitMean = inla.mesh.project(mapBasis, 
                              model_redfish$summary.random$field$mean[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  fit.025 = inla.mesh.project(mapBasis, 
                              model_redfish$summary.random$field$`0.025quant`[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  fit.975 = inla.mesh.project(mapBasis, 
                              model_redfish$summary.random$field$`0.975quant`[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  
  # Build maps with confidence intervals
  mapMean[[i]] = raster(t(fitMean[,ncol(fitMean):1]),
                         xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                         ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                         crs = crs(gulf))
  
  map.025[[i]] = raster(t(fit.025[,ncol(fit.025):1]),
                         xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                         ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                         crs = crs(gulf))
  
  map.975[[i]] = raster(t(fit.975[,ncol(fit.975):1]),
                         xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                         ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                         crs = crs(gulf))
  
}

# Convert to a rasterStack
rasterMean = stack(mapMean)
raster.025 = stack(map.025)
raster.975 = stack(map.975)


for(i in 1:5){
  values(rasterMean)[,i] = as.vector(t(mapMean[[i]]))
  values(raster.025)[,i] = as.vector(t(map.025[[i]]))
  values(raster.975)[,i] = as.vector(t(map.975[[i]]))
}

# Mask the region
rasterMeanMask = mask(rasterMean, gulf)
raster.025Mask = mask(raster.025, gulf)
raster.975Mask = mask(raster.975, gulf)

# Plot the results
par(mfrow = c(5,3), mar = c(1,1,3,1))

zlimRange = range(values(raster.025Mask[[i]]),
                   values(rasterMeanMask[[i]]),
                   values(raster.975Mask[[i]]), na.rm = TRUE)

for(i in 1:5){
  plot(raster.025Mask[[i]],
       zlim = zlimRange,
       main = "2.5%")
  
  plot(rasterMeanMask[[i]],
       zlim = zlimRange,
       main = "Mean")
  
  plot(raster.975Mask[[i]],
       zlim = zlimRange,
       main = "97.5%")
  
}
