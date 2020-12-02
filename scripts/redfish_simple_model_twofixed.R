whitehake ~ 0 + bTemp + bSal +
  f(field, # temporal autocorr
    model=SPDE, # spatial autocorr structure
    group = field.group, # how they're grouped
    control.group=list(model='ar1', hyper=hSpec) # build with model parameters
  ) +
  f(
    strat,
    model = "iid",
    hyper = list(theta = list(prior = "gaussian",
                              param = c(0, 0.001)))
  )


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
                         cutoff=maxEdge,
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
Alist_simple = as.list(rep(1,4))
Alist_simple[[1]] = A

effect_simple = list(Field,
                     bTemp =  explan@data$bottom.temperature,
                     bSal = explan@data$bottom.salinity,
                     stratum = explan@data$stratum)

Stack_redfish_simple = inla.stack(data=list(redfish = number@data$`REDFISH UNSEPARATED`),
                                  A = Alist_simple,
                                  effects = effect_simple,
                                  tag="basis")

form_redfish_simple_1 = redfish ~ 0 + bTemp + bSal +
  f(field, # temporal autocorr
    model=SPDE, # spatial autocorr structure
    group = field.group, # how they're grouped
    control.group=list(model='ar1', hyper=hSpec) # build with model parameters
  ) + f(stratum,
    model = "iid",
    hyper = list(theta = list(prior = "gaussian",
                              param = c(0, 0.001)))
  )


model_redfish_simple_rand_rand_1 = inla(form_redfish_simple_1,
                                      data = inla.stack.data(Stack_redfish_simple),
                                      family="nbinomial",
                                      control.family =list(link="log", hyper = list(theta=precPrior)),
                                      control.predictor=list(A=inla.stack.A(Stack_redfish_simple),
                                                             compute=TRUE, link = 1),
                                      control.compute=list(waic=TRUE),
                                      verbose = FALSE)

summary(model_redfish_simple_rand_rand_1)
saveRDS(model_redfish_simple_rand_rand_1, 
               here('./model-output/model_redfish_simple_rand_space-stratum.rds'))

ExpY <- model_redfish_simple_rand_rand$summary.fitted.values[1:nrow(number@data),"mean"]

#For the variance we also need the posterior mean value of θ, which is obtained via
phi1 <- model_redfish_simple_rand_rand$summary.hyper[1, "mean"]

#The variance is then given by
VarY <- ExpY * (1 - ExpY) / (1 + phi1)

# And once we have the mean and the variance we can calculate Pearson residuals.
E1 <- (number@data$`REDFISH UNSEPARATED` - ExpY) / sqrt(VarY)

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = model_redfish_simple_rand_rand$summary.fitted.values$mean[1:nrow(number@data)], 
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     xlim = c(0, 1))
abline(h = 0, lty = 2)


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
                              model_redfish_simple_rand$summary.random$field$mean[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  fit.025 = inla.mesh.project(mapBasis, 
                              model_redfish_simple_rand$summary.random$field$`0.025quant`[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  fit.975 = inla.mesh.project(mapBasis, 
                              model_redfish_simple_rand$summary.random$field$`0.975quant`[1:SPDE$n.spde + (i - 1) * SPDE$n.spde])
  
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
