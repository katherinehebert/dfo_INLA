# Load R package
library(INLA)
library(raster)
library(animation)
library(classInt)
library(rgeos)

#___________
# Load model 
#___________
model <- readRDS("Snailfish unidentified - model.RDS")
ID <- readRDS("Snailfish unidentified - ID.RDS") # Important for projection

#model <- readRDS("Atlantic cod - model.RDS")
#ID <- readRDS("Atlantic cod - ID.RDS") # Important for projection

#__________________
# Load mapping info
#__________________
gulf <- readRDS("./Data/gulfSampled.RDS")
canada <- readRDS("./Data/Canada.RDS")

# Load data (to draw strata)
explan <- readRDS("./Data/explan.RDS")

# Independent year
year <- 1971:2020
nyear <- length(year)

#_____________________________
# Load mesh and related object
#_____________________________
meshSpace <- readRDS("./Mesh/meshSpace.RDS")

############
# Plot model
############
# Build basis of the map
stepsize <- 1000
rangeX <- range(meshSpace$loc[,1])
rangeY <- range(meshSpace$loc[,2])
nxy <- round(c(diff(rangeX), 
               diff(rangeY)) / stepsize)

# Define basis of the map (warning is OK)
mapBasis <- inla.mesh.projector(meshSpace,
                                xlim = rangeX,
                                ylim = rangeY,
                                dim = nxy,
                                crs = meshSpace$crs)

# Result objects
mapGulfMean <- vector("list", length = nyear)
mapGulf.025 <- vector("list", length = nyear)
mapGulf.975 <- vector("list", length = nyear)

for(j in 1:nyear){
  if(!is.null(model[[j]])){
    # Model prediction
    fitMean <- inla.mesh.project(mapBasis, 
                                 model[[j]]$summary.fitted.values$mean[ID[,j]])
    fit.025 <- inla.mesh.project(mapBasis, 
                                 model[[j]]$summary.fitted.values$`0.025quant`[ID[,j]])
    fit.975 <- inla.mesh.project(mapBasis, 
                                 model[[j]]$summary.fitted.values$`0.975quant`[ID[,j]])
  }else{
    # Model prediction
    fitMean <- inla.mesh.project(mapBasis, 
                                 rep(0, meshSpace$n))
    fit.025 <- inla.mesh.project(mapBasis, 
                                 rep(0, meshSpace$n))
    fit.975 <- inla.mesh.project(mapBasis, 
                                 rep(0, meshSpace$n))
  }  
  # Build maps with confidence intervals
  mapMean <- raster(t(fitMean[,ncol(fitMean):1]),
                    xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                    ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                    crs = crs(gulf))
  
  map.025 <- raster(t(fit.025[,ncol(fit.025):1]),
                    xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                    ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                    crs = crs(gulf))
  
  map.975 <- raster(t(fit.975[,ncol(fit.975):1]),
                    xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                    ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                    crs = crs(gulf))
  
  # Crop and mask
  mapGulfMean[[j]] <- mask(crop(mapMean, extent(gulf)), gulf)
  mapGulf.025[[j]] <- mask(crop(map.025, extent(gulf)), gulf)
  mapGulf.975[[j]] <- mask(crop(map.975, extent(gulf)), gulf)
}

# Make stack
mapGulfStackMean <- stack(mapGulfMean)
mapGulfStack.025 <- stack(mapGulf.025)
mapGulfStack.975 <- stack(mapGulf.975)

# Change layers name
names(mapGulfStackMean) <- paste0("y", year)
names(mapGulfStack.025) <- paste0("y", year)
names(mapGulfStack.975) <- paste0("y", year)

# Extract most extreme values 
val <- c(values(mapGulfStack.025),
         values(mapGulfStack.975))

# Define class intervals
fixBreak <- c(0,exp(seq(0,9.25, length = 99)), max(val, na.rm = TRUE))

# Warning OK
breaks <- classIntervals(val,
                         style = "fixed",
                         fixedBreaks = fixBreak)

color <- colorRampPalette(c("grey90", "steelblue4", "steelblue2", 
                            "steelblue1", "gold", "red1", "red4"),
                          bias = 1)


# Plot figure for model mean
saveGIF(
  {
    for(j in 1:nyear){
      # Base of the figure
      par(mfrow = c(1,4), mar=c(1,1,1,1))
      #=====
      # 2.5%
      #=====
      # Base of plot
      plot(mapGulfStack.025@extent,
           xlim = c(mapGulfStack.025@extent@xmin,
                    mapGulfStack.025@extent@xmax),
           ylim = c(mapGulfStack.025@extent@ymin,
                    mapGulfStack.025@extent@ymax),
           type = "n", axes=FALSE,
           xlab = "", ylab = "", asp = 1)
      
      mtext(text = year[j],side = 3,outer = TRUE, cex = 3)
      
      # Background
      rect(mapGulfStack.025@extent@xmin-25000,
           mapGulfStack.025@extent@ymin-25000,
           mapGulfStack.025@extent@xmax+25000,
           mapGulfStack.025@extent@ymax+25000,
           col = rgb(193/256,236/256,250/256),
           border = rgb(193/256,236/256,250/256))
      
      # Model
      plot(mapGulfStack.025[[j]], 
           col = color(length(breaks$brks)),
           add = TRUE,
           breaks = breaks$brks,
           legend = FALSE)
      
      plot(canada, add = TRUE, col ="darkkhaki")
      legend("topright", "2.5%", bty = "n", cex = 3, adj = c(0,0))
      
      #=====
      # Mean
      #=====
      # Base of plot
      plot(mapGulfStackMean@extent,
           xlim = c(mapGulfStackMean@extent@xmin,
                    mapGulfStackMean@extent@xmax),
           ylim = c(mapGulfStackMean@extent@ymin,
                    mapGulfStackMean@extent@ymax),
           type = "n", axes=FALSE,
           xlab = "", ylab = "", asp = 1)
      
      # Background
      rect(mapGulfStackMean@extent@xmin-25000,
           mapGulfStackMean@extent@ymin-25000,
           mapGulfStackMean@extent@xmax+25000,
           mapGulfStackMean@extent@ymax+25000,
           col = rgb(193/256,236/256,250/256),
           border = rgb(193/256,236/256,250/256))
      
      # Model
      plot(mapGulfStackMean[[j]], 
           col = color(length(breaks$brks)),
           add = TRUE,
           breaks = breaks$brks,
           legend = FALSE)

      plot(canada, add = TRUE, col ="darkkhaki")
      legend("topright", "Mean", bty = "n", cex = 3, adj = c(0,0))
      
      #======
      # 97.5%
      #======
      # Base of plot
      plot(mapGulfStack.975@extent,
           xlim = c(mapGulfStack.975@extent@xmin,
                    mapGulfStack.975@extent@xmax),
           ylim = c(mapGulfStack.975@extent@ymin,
                    mapGulfStack.975@extent@ymax),
           type = "n", axes=FALSE,
           xlab = "", ylab = "", asp = 1)
      
      # Background
      rect(mapGulfStack.975@extent@xmin-25000,
           mapGulfStack.975@extent@ymin-25000,
           mapGulfStack.975@extent@xmax+25000,
           mapGulfStack.975@extent@ymax+25000,
           col = rgb(193/256,236/256,250/256),
           border = rgb(193/256,236/256,250/256))
      # Model
      plot(mapGulfStack.975[[j]], 
           col = color(length(breaks$brks)),
           add = TRUE,
           breaks = breaks$brks,
           legend = FALSE)

      plot(canada, add = TRUE, col ="darkkhaki")
      legend("topright", "97.5%", bty = "n", cex = 3, adj = c(0,0))
      
      # Legend
      par(mar =c(5,28,10,28))
      image(matrix(1:length(fixBreak), nrow = 1), col = color(length(fixBreak)),
            xaxt = "n", yaxt = "n")
      axis(4, at = seq(0.05,0.95,by=0.1),
           #           labels = round(fixBreak[seq(5,95,by=10)]), 
           labels = c(1, 3, 10, 25, 60, 150, 400, 1000, 2500, 6500), 
           las = 1, cex.axis = 2.5)
      title(main = year[j], cex.main = 5)
      
    }
  }, 
  movie.name = paste0("Snailfish unidentified.gif"),
  ani.height = 500,
  ani.width = 2400
)
