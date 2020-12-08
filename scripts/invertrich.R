###########################
###########################
## This script will build, run, and summarise **simple** hypothesis models in 
## INLA for invert richness for the DFO St. Lawrence dataset
###########################
###########################
## date: 2020-12-08
## author: Katherine Hébert
###########################
###########################

# The hypotheses models are:

# 1) btemp (energy hypothesis)
# 2) ice cover, bSalinity (productivity hypothesis)
# 3) bTemp + sTemp (climate stability) # note... this is exactly the same as energy hypothesis if we remove sTemp
# 4) slope, bSalinity (spatial heterogeneity hypothesis)
# 5) depth (stress hypothesis)

# full model: invertrich ~ bTemp + sTemp + ice cover + bSalinity + slope + depth

# set-up =======================================================================

library(INLA)
library(tidyverse)
library(ggregplot) # devtools::install_github("gfalbery/ggregplot")

# load data
explan <- readRDS("data/explan_ice.RDS")
invert <- readRDS("data/invertRich.RDS")


# prepare data for model =======================================================

# build A matrix that will hold parameters to construct the model
Alist <- as.list(rep(1,6)) # list with single values

# organise effects 
effect <- list(bTemp =  explan@data$bottom.temperature,
               sTemp = explan@data$surface.temperature,
               bSal = explan@data$bottom.salinity,
               slope = explan@data$slope,
               bDep = explan@data$depth,
               ice = explan@data$ice.duration
)

# build stack
Stack <- inla.stack(data=list(invert = invert@data$species.invertebrate.number),
                    A = Alist,
                    effects = effect,
                    tag="basis")

# set prior
precPrior <- list(prior='pc.prec', param=c(1, 0.01))


# build models =================================================================

# function to use different forms but keep all other specifications the same
quick_inla <- function(form){
  inla(form, 
       data = inla.stack.data(Stack),
       family="nbinomial",
       control.family = list(link="log", hyper = list(theta=precPrior)),
       control.predictor = list(A = inla.stack.A(Stack), 
                                compute = TRUE, 
                                link = 1),
       control.compute = list(waic = TRUE, config = TRUE))
}

# full model ----
full <- invert ~ 0 + bTemp + bSal + slope + bDep + ice
m_full <- quick_inla(full)

# energy hypothesis ----
energy <- invert ~ 0 + bTemp 
m_energy <- quick_inla(energy)

# productivity hypothesis ----
prod <- invert ~ 0 + bSal + ice
m_prod <- quick_inla(prod)

# climate stability hypothesis ----
clim <- invert ~ 0 + bTemp + sTemp
m_clim <- quick_inla(clim)

# spatial heterogeneity hypothesis ----
sphet <- invert ~ 0 + slope + bSal
m_sphet <- quick_inla(sphet)

# stress hypothesis ----
stress <- invert ~ 0 + bDep
m_stress <- quick_inla(stress)

# save models
saveRDS(m_full, "model-output/invertrich_full.rds")
saveRDS(m_energy, "model-output/invertrich_energy.rds")
saveRDS(m_prod, "model-output/invertrich_prod.rds")
saveRDS(m_clim, "model-output/invertrich_clim.rds")
saveRDS(m_sphet, "model-output/invertrich_sphet.rds")
saveRDS(m_stress, "model-output/invertrich_stress.rds")


# model summary ================================================================

Efxplot(ModelList = list(m_full, m_energy, m_prod, m_clim, m_sphet, m_stress),
        ModelNames = c("Full", "Energy", "Productivity", "Climate stability", "Spatial heterogeneity", "Stress"),
        Size = 2) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 11))
ggsave("figures/invertrich_coefplot.pdf", width = 5.59, height = 4.84)

# model selection ==============================================================

# extract waic to compare
waic_df <- data.frame(model = c("full", "energy", "prod", "clim", "sphet", "stress"), 
                      waic = NA)
waic_df$waic <- c(m_full$waic$waic,
                  m_energy$waic$waic,
                  m_prod$waic$waic,
                  m_clim$waic$waic,
                  m_sphet$waic$waic,
                  m_stress$waic$waic
)
waic_df <- waic_df[order(waic_df$waic),]
waic_df$model <- factor(waic_df$model, levels = waic_df$model)

# plot
ggplot(waic_df) + 
  geom_point(aes(x = model, y = waic), size  = 3) + 
  theme_minimal() + theme(axis.title.y = element_text(size = 14),
                          axis.text.x = element_text(size = 12)) +
  labs(x = "", title = "Model comparison")

# save outputs
ggsave("figures/invertrich_waic.pdf", width = 5.57, height = 3.37)
saveRDS(waic_df, "outputs/invertrich_waic.rds")