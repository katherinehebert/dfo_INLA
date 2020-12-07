# Script to get ice cover data

# https://github.com/duplisea/gslea
# devtools::install_github("duplisea/gslea", build_vignettes = TRUE)
library(gslea)

# EAR = Regions 5 and 6 covers the entire sGSL

# Duration of the ice season (number of days)
EA.plot.f(variables = find.vars.f("ice.duration"), years = 1800:2100, EARs=5:6)
ice_duraction <- EA.query.f(find.vars.f("ice.duration"), years = 1800:2100, EARs=5:6)

# to do:
# match with explanatory variables