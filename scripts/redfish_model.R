####practice script for ciee work
library(here)
library(raster)
library(INLA)
library(animation)

fishRich = readRDS(here('fishRich.RDS'))
catchTotal = readRDS(here('catchTotal.RDS'))
explan = readRDS(here('explan.RDS'))
gulf = readRDS(here('gulf.RDS'))
number = readRDS(here('number.RDS'))

number@data$`WHITE HAKE`

redfish = number@data$`REDFISH UNSEPARATED`

maxEdge = 7500
meshSpace = inla.mesh.2d(boundary = gulf, 
                         max.edge = maxEdge*c(0.5,2),
                         offset = c(100000,20000),
                         cutoff = maxEdge/2)
par(mar = c(1,1,1,1))
plot(meshSpace, main = "", asp = 1)

meshSpace$n

time = as.numeric(fishRich@endTime)/60
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
             verbose = TRUE)

summary(model_redfish)