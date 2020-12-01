# Load R package
library(spacetime)
library(INLA)

#__________
# Load data
#__________
explan <- readRDS("explan.RDS")
gulf <- readRDS("gulf.RDS")

#_____________________
## Reorganise the data
#_____________________
### Space
xy <- coordinates(explan)
colnames(xy) <- c("locx", "locy")

### Time
time <- as.numeric(explan@endTime)/60 # Number of minutes since Jan 1, 1970 
timeBasis <- time - 884844 # Time starting at 1

#______________
## Build meshes
#______________
### Temporal mesh
timeGr <- 11
timeSeq <- seq(min(timeBasis),
               max(timeBasis),
               length = timeGr)

meshTime <- inla.mesh.1d(timeSeq)

### Spatial mesh
maxEdge <- 35000
meshSpace <- inla.mesh.2d(boundary = gulf, 
                          max.edge=maxEdge * c(0.5,2),
                          cutoff=maxEdge/4,
                          offset = c(10000,20000))

plot(meshSpace, asp = 1)
