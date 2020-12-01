# Load R package
library(spacetime)
library(INLA)

#__________
# Load data
#__________
explan <- readRDS("explan.RDS")
number <- readRDS("number.RDS")

#__________________________________________
# Model with a random effect on the stratum
#__________________________________________
model <- inla(CAPELIN ~ f(explan@data$stratum,
                          model = "iid",
                          hyper = list(theta = list(prior = "gaussian",
                                                    param = c(0, 0.001)))),
              data = number@data, 
              family = "nbinomial")

# Extract result
summary(model)
model$summary.hyperpar
