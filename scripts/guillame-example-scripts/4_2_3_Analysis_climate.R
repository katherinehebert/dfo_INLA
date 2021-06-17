# Load R package
library(spacetime)
library(INLA)
library(raster)

#__________
# Load data
#__________
# Data
sp <- readRDS("./Data/number.RDS")
explan <- readRDS("./Data/explan.RDS")
explanMesh <- readRDS("./Data/explanMesh.RDS")

# Focus only on species with 50 occurrences or more across full data
sp@data <- sp@data[,which(colSums(sp@data) >= 50)]

# WAIC results from hypothesis
WAIC <- read.csv("./Results/WAICDiff_Hypothese.csv", row.names = 1)

#__________________
# Load mapping info
#__________________
gulf <- readRDS("./Data/gulf.RDS")
canada <- readRDS("./Data/Canada.RDS")

#_____________________________
# Load mesh and related object
#_____________________________
meshSpace <- readRDS("./Mesh/meshSpace.RDS")
barrierTriangle <- readRDS("./Mesh/barrierTriangle.RDS")

#_________________________________________________
# Space only barrier model across years separately
#_________________________________________________
# Get species info to build model
year <- unique(as.POSIXlt(sp@endTime)$year)
nyear <- length(year)

# Species to model
i = 64

# result object
model <- vector("list", length = nyear)

# ID for prediction
ID <- matrix(NA,nrow = meshSpace$n, ncol = nyear)

for(j in 1:nyear){
  pointerYear <-  which(as.POSIXlt(sp@endTime)$year == year[j])
  locYear <- coordinates(sp@sp)[pointerYear,]
  spYear <- sp@data[pointerYear,i]
  sum(spYear)

  if(sum(spYear)>0){

    pointerExplanMesh <- which(as.POSIXlt(explanMesh@endTime)$year == year[j])

    # For fish species
    if(WAIC$fishInver[i] == "fish"){
      # Build stack for estimation
      AEst <- inla.spde.make.A(meshSpace, locYear)
      stackEst <- inla.stack(data = list(y = spYear),
                             A = list(AEst, 1,1),
                             effect = list(s = 1:meshSpace$n,
                                           intercept = rep(1, length(pointerYear)),
                                           surface.temperature.5.var = scale(explan@data$surface.temperature.5.var)[pointerYear]),
                             tag = "est")

      # Build stack for prediction
      APred <-  inla.spde.make.A(meshSpace, meshSpace$loc[,1:2])
      stackPred <- inla.stack(data = list(y = NA),
                              A = list(APred,1,1),
                              effect = list(s = 1:meshSpace$n,
                                            intercept = rep(1, length(pointerExplanMesh)),
                                            surface.temperature.5.var = scale(explanMesh@data$surface.temperature.5.var)[pointerExplanMesh]),
                              tag = "pred")

      # Formula
      Formula <- y ~ 0 + intercept + surface.temperature.5.var + f(s, model = barrierModel)
    }

    # For invertebrate species
    if(WAIC$fishInver[i] == "invertebrate"){
      # Build stack for estimation
      AEst <- inla.spde.make.A(meshSpace, locYear)
      stackEst <- inla.stack(data = list(y = spYear),
                             A = list(AEst, 1,1),
                             effect = list(s = 1:meshSpace$n,
                                           intercept = rep(1, length(pointerYear)),
                                           bottom.temperature.5.var = scale(explan@data$bottom.temperature.5.var)[pointerYear]),
                             tag = "est")

      # Build stack for prediction
      APred <-  inla.spde.make.A(meshSpace, meshSpace$loc[,1:2])
      stackPred <- inla.stack(data = list(y = NA),
                              A = list(APred,1,1,1),
                              effect = list(s = 1:meshSpace$n,
                                            intercept = rep(1, length(pointerExplanMesh)),
                                            bottom.temperature.5.var = scale(explanMesh@data$bottom.temperature.5.var)[pointerExplanMesh]),
                              tag = "pred")

      # Formula
      Formula <- y ~ 0 + intercept + bottom.temperature.5.var + f(s, model = barrierModel)
    }

    # Join stacks
    stack <- inla.stack(stackEst, stackPred)

    # ID in stack for prediction
    ID[,j] <- inla.stack.index(stack, tag="pred")$data

    # Build basis of barrier model
    barrierModel <- inla.barrier.pcmatern(meshSpace,
                                          barrier.triangles = barrierTriangle,
                                          prior.range = c(100000, 0.5),
                                          prior.sigma = c(30, 0.5))

    # PC prior for standard deviation
    sdPCprior <- list(prior = "pc.prec", param = c(20, 0.5))
    # * Need to define prior carefully

    # Fit model
    model[[j]] <- inla(Formula,
                       data = inla.stack.data(stack),
                       control.predictor = list(A = inla.stack.A(stack),
                                                compute = TRUE,
                                                link = 1),
                       family = "nbinomial",
                       control.inla = list(int.strategy = "eb"),
                       control.family = list(hyper = list(theta = sdPCprior)),
                       control.compute = list(waic = TRUE))
  }
  print(paste(1900 + year[j],"Done"))
}

saveRDS(model, file = paste0(colnames(sp@data)[i]," - model.RDS"))
saveRDS(ID, file = paste0(colnames(sp@data)[i]," - ID.RDS"))
