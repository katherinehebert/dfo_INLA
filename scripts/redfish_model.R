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
library(patchwork)
#library(ggregplot)


#read in all data
source(here('./scripts/read_in_all_data.R'))

### Space
xy = coordinates(explan)
colnames(xy) = c("locx", "locy")

### Time
time = as.numeric(explan@endTime)/60 # Number of minutes since Jan 1, 1970 
timeBasis = time - 884844 # Time starting at 1

#______________
## Build meshes
#______________
### Temporal mesh
timeGr = 11
timeSeq = seq(min(timeBasis),
               max(timeBasis),
               length = timeGr)

meshTime = inla.mesh.1d(timeSeq)

### Spatial mesh
maxEdge = 35000
meshSpace = inla.mesh.2d(boundary = gulf, 
                          max.edge=maxEdge * c(0.5,2),
                          cutoff=maxEdge/4,
                          offset = c(10000,20000))

plot(meshSpace, asp = 1)

SPDE = inla.spde2.pcmatern(mesh=meshSpace,
                           alpha=2,
                           prior.range=c(1000, 0.5),
                           prior.sigma=c(35, 0.5))


hSpec = list(theta=list(prior='pccor1', param=c(0.2, 0.9)))
precPrior = list(prior='pc.prec', param=c(1, 0.01)) 

Field = inla.spde.make.index("field", n.spde=SPDE$n.spde,
                             n.group = meshTime$n)


# For estimation
A = inla.spde.make.A(meshSpace,
                     loc=coordinates(number@sp),
                     group = timeBasis,
                     group.mesh = meshTime)

######## model fitting for some other models to compare

#include depth & random effect for the full model
Alist_btbsdp_rand = as.list(rep(1,5))
Alist_btbsdp_rand[[1]] = A

effect_btbsdp_rand = list(Field,
                   bTemp =  explan@data$bottom.temperature,
                   bSal = explan@data$bottom.salinity,
                   depth = explan@data$depth,
                   stratum = explan@data$stratum)

Stack_redfish_btbdp_rand = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                               A = Alist_btbsdp_rand,
                               effects = effect_btbsdp_rand,
                               tag="basis")

form_redfish_btbdp_rand = redfish ~ 0 + bTemp + bSal + depth + f(explan@data$stratum,
                                                                 model = "iid",
                                                                 hyper = list(theta = list(prior = "gaussian",
                                                                                           param = c(0, 0.001))))
  

model_redfish_btbdp_rand = inla(form_redfish_btbdp_rand,
                         data = inla.stack.data(Stack_redfish_btbdp_rand),
                         family="nbinomial",
                         control.family =list(link="identity", hyper = list(theta=precPrior)),
                         control.predictor=list(A=inla.stack.A(Stack_redfish_btbdp_rand),
                                                compute=TRUE, link = 1),
                         control.compute=list(waic=TRUE),
                         verbose = FALSE)

summary(model_redfish_btbdp_rand)

# model with all fixed but no random
Alist_btbsdp = as.list(rep(1,5))
Alist_btbsdp[[1]] = A

effect_btbsdp = list(Field,
                          bTemp =  explan@data$bottom.temperature,
                          bSal = explan@data$bottom.salinity,
                          depth = explan@data$depth,
                     intercept = rep(1,nrow(explan@data)))

Stack_redfish_btbdp = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                                      A = Alist_btbsdp,
                                      effects = effect_btbsdp,
                                      tag="basis")

form_redfish_btbdp = redfish ~ 0 + bTemp + bSal + depth + intercept +
  f(field, model=SPDE,
    group = field.group,
    control.group=list(model='ar1', hyper=hSpec))

model_redfish_btbdp = inla(form_redfish_btbdp,
                                data = inla.stack.data(Stack_redfish_btbdp),
                                family="gaussian",
                                control.family =list(link="identity", hyper = list(theta=precPrior)),
                                control.predictor=list(A=inla.stack.A(Stack_redfish_btbdp),
                                                       compute=TRUE, link = 1),
                                control.compute=list(waic=TRUE),
                                verbose = FALSE)

summary(model_redfish_btbdp)

ExpY <- model_redfish_btbdp$summary.fitted.values[1:nrow(number@data),"mean"]

#For the variance we also need the posterior mean value of θ, which is obtained via
phi1 <- model_redfish_btbdp$summary.hyper[1, "mean"]

#The variance is then given by
VarY <- ExpY * (1 - ExpY) / (1 + phi1)

# And once we have the mean and the variance we can calculate Pearson residuals.
E1 <- (number@data$`REDFISH UNSEPARATED` - ExpY) / sqrt(VarY)

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = model_redfish_btbdp$summary.fitted.values$mean[1:nrow(number@data)], 
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     xlim = c(0, 1))
abline(h = 0, lty = 2)
# model fitting for the basic model we did in the group (no random no depth)
Alist_btbs = as.list(rep(1,3))
Alist_btbs[[1]] = A

effect_btbs = list(Field,
                   bTemp =  explan@data$bottom.temperature,
                   bSal = explan@data$bottom.salinity)

Stack_redfish_btbs = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                               A = Alist_btbs,
                               effects = effect_btbs,
                               tag="basis")


form_redfish_btbs = redfish ~ 0 + bTemp + bSal +
  f(field, model=SPDE,
    group = field.group,
    control.group=list(model='ar1', hyper=hSpec))

model_redfish_btbs = inla(form_redfish_btbs,
                         data = inla.stack.data(Stack_redfish_btbs),
                         family="gaussian",
                         control.family =list(link="identity", hyper = list(theta=precPrior)),
                         control.predictor=list(A=inla.stack.A(Stack_redfish_btbs),
                                                compute=TRUE, link = 1),
                         control.compute=list(waic=TRUE),
                         verbose = FALSE)

summary(model_redfish_btbs)

# model with bottom temperature and depth 
Alist_btdp = as.list(rep(1,3))
Alist_btdp[[1]] = A

effect_btdp = list(Field,
                   bTemp =  explan@data$bottom.temperature,
                   depth = explan@data$depth)

Stack_redfish_btdp = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                               A = Alist_btdp,
                               effects = effect_btdp,
                               tag="basis")


form_redfish_btdp = redfish ~ 0 + bTemp + depth +
  f(field, model=SPDE,
    group = field.group,
    control.group=list(model='ar1', hyper=hSpec))

model_redfish_btdp = inla(form_redfish_btdp,
                         data = inla.stack.data(Stack_redfish_btdp),
                         family="gaussian",
                         control.family =list(link="identity", hyper = list(theta=precPrior)),
                         control.predictor=list(A=inla.stack.A(Stack_redfish_btdp),
                                                compute=TRUE, link = 1),
                         control.compute=list(waic=TRUE),
                         verbose = FALSE)

summary(model_redfish_btdp)

# model with depth and bottom salinity
Alist_bsdp = as.list(rep(1,3))
Alist_bsdp[[1]] = A

effect_bsdp = list(Field,
                   bSal =  explan@data$bottom.salinity,
                   depth = explan@data$depth)

Stack_redfish_bsdp = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                                A = Alist_bsdp,
                                effects = effect_bsdp,
                                tag="basis")


form_redfish_bsdp = redfish ~ 0 + bTemp + depth +
  f(field, model=SPDE,
    group = field.group,
    control.group=list(model='ar1', hyper=hSpec))

model_redfish_bsdp = inla(form_redfish_bsdp,
                          data = inla.stack.data(Stack_redfish_bsdp),
                          family="gaussian",
                          control.family =list(link="identity", hyper = list(theta=precPrior)),
                          control.predictor=list(A=inla.stack.A(Stack_redfish_bsdp),
                                                 compute=TRUE, link = 1),
                          control.compute=list(waic=TRUE),
                          verbose = FALSE)

summary(model_redfish_bsdp)



#Extracting fitted values and residuals
#Once a model has been executed model validation has to be applied and
#for this we need fitted values, variances and Pearson residuals. To get the
#fitted values and Pearson residuals we need to follow the expression in
#Equation (23.1) or (23.3). The expected values μ
#i are called ExpY in the
#code below and are obtained via
#ExpY <- I2a$summary.fitted.values[1:N,"mean"]
#For the variance we also need the θ, which is obtained via
#phi1 <- I2a$summary.hyper[1, "mean"]
#The posterior mean value of θ is 12.86. The variance is then given by
#VarY <- ExpY * (1 - ExpY) / (1 + phi1)
#And once we have the mean and the variance we can calculate Pearson
#residuals.
#E1 <- (CR2$CCAPropTran - ExpY) / sqrt(VarY)

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
                              model_redfish_btbdp$summary.random$field$mean[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  fit.025 = inla.mesh.project(mapBasis, 
                              model_redfish_btbdp$summary.random$field$`0.025quant`[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  fit.975 = inla.mesh.project(mapBasis, 
                              model_redfish_btbdp$summary.random$field$`0.975quant`[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  
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

# convert to data frames for ggplot
raster.025Mask_df <- as.data.frame(raster.025Mask, xy = TRUE, na.rm = TRUE)
rasterMeanMask_df <- as.data.frame(rasterMeanMask, xy = TRUE, na.rm = TRUE)
raster.975Mask_df <- as.data.frame(raster.975Mask, xy = TRUE, na.rm = TRUE)


# function to plot model results using ggplot2
result_map <- function(raster_df, column_name, title){
  p <- ggplot() +
    geom_raster(data = raster_df, 
                aes(x = x, y = y, fill = get(column_name))) + 
    labs(fill = "Abundance") +
    scale_fill_distiller(palette = "RdBu", 
                         limits = c(-max(abs(raster.975Mask_df[,column_name])), 
                                    max(abs(raster.975Mask_df[,column_name])))) +
    theme_void()  
  return(p)
}

# 2.5% credible limit
p1.025 <- result_map(raster.025Mask_df, "layer.1") + theme(legend.position = "none")
p2.025 <- result_map(raster.025Mask_df, "layer.2") + theme(legend.position = "none")

# mean
p1.mean <- result_map(rasterMeanMask_df, "layer.1") + theme(legend.position = "none")
p2.mean <- result_map(rasterMeanMask_df, "layer.2") + theme(legend.position = "none")

# 97.5% credible limit
p1.975 <- result_map(raster.975Mask_df, "layer.1") + theme(legend.position = "bottom")
p2.975 <- result_map(raster.975Mask_df, "layer.2") + theme(legend.position = "bottom")

# patchwork together
(p1.025 + p2.025) / (p1.mean + p2.mean) / (p1.975 + p2.975)





###################################### how to do random effects:
form_redfish_btbdp = redfish ~ 0 + bTemp + bSal + depth + intercept +
  f(field, model=SPDE,
    group = field.group,
    control.group=list(model='ar1', hyper=hSpec))
#each different level of our random sample is a random effect - similar to mgcv
form = y ~ f(ID = 1:nrow(explan@data), model = "iid")

model_redfish_btbdp = inla(form_redfish_btbdp,
                           data = inla.stack.data(Stack_redfish_btbdp),
                           family="gaussian",
                           control.family =list(link="identity", hyper = list(theta=precPrior)),
                           control.predictor=list(A=inla.stack.A(Stack_redfish_btbdp),
                                                  compute=TRUE, link = 1),
                           control.compute=list(waic=TRUE),
                           verbose = FALSE)
