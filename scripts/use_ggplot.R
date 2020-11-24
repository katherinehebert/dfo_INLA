# Script to explore the data using sf and ggplot

library(sf)
library(tidyverse)
require(here)

# load data
source(here("./scripts/read_in_all_data.R"))

# function to add coordinates to dataset in preparation for conversion to sf
add_coords <- function(stidf_object){
  # assign datasets to objects to simplify the code a bit
  df <- stidf_object@data 
  coords <- stidf_object@sp@coords

  # add coordinate columns to dataset
  df <- mutate(df,
           longitude = coords[,1],
           latitude = coords[,2],
           time = lubridate::date(stidf_object@endTime),
           year = lubridate::year(stidf_object@endTime)
           )
  return(df)
}

# convert datasets to sf objects ----

number_sf <- add_coords(number) %>%
  st_as_sf(x = ., coords = c("longitude", "latitude"))
st_crs(number_sf) <- number@sp@proj4string # set coordinate ref system

explan_sf <- add_coords(explan) %>%
  st_as_sf(x = ., coords = c("longitude", "latitude"))
st_crs(explan_sf) <- explan@sp@proj4string # set coordinate ref system

gulf_sf <- st_as_sf(gulf)
st_crs(gulf_sf) <- gulf@proj4string # set coordinate ref system

# maps -----

# ggplot map of explanatory variable
explanatory_variable_map = function(variable, years) {
  
  # Function uses predefined data of the gulf and the explan_sf df,but takes in
  # the variable of interest to plot and returns a map
  #
  # Parameters:
  #            variable (TYPE = character): vector (dataframe column) of interest 
  #                                      passed as the name of the column
  #            years (TYPE = char vector): concatenated vector of the study 
  #                                       years of interested passed as integers
  # Returns:
  #            plot (TYPE = ggplot object): plot of desired variable
  #
  # EXAMPLE USAGE:
  # > explanatory_variable_map('bottom.temperature', c(1971, 1981, 1991))
  #
  
  return(
    
    ggplot() +
      geom_sf(data = gulf_sf, fill = "white") +
      geom_sf(data = filter(explan_sf, 
                            # pick a sequence of years to check out
                            year %in% years),
              aes_string(col = variable)) +
      facet_wrap(~year) +
      scale_color_viridis_c(option = "viridis") +
      labs(color = gsub("\\.", "\n", variable)) +
      theme_void() 
  )
  
}

explanatory_variable_map('bottom.temperature', c(1971, 1981, 1991, 
                                                 2001, 2011, 2019))
explanatory_variable_map('surface.temperature', c(1971, 1981, 1991, 
                                                 2001, 2011, 2019))


## ggplot map of abundance for a given species
species_abundance_map = function(species, years) {
  
  # Function uses predefined data of the gulf and the number_sf df,but takes in
  # the species of interest to plot and returns a map of it's abundance
  #
  # Parameters:
  #            species (TYPE = character): species name of interest 
  #                                      passed as the name of the column
  #            years (TYPE = char vector): concatenated vector of the study 
  #                                       years of interested passed as integers
  # Returns:
  #            plot (TYPE = ggplot object): plot of desired variable
  #
  # EXAMPLE USAGE:
  # > species_abundance_map('WHITE HAKE', c(1971, 1981, 1991))
  #
  
  return(
   
    ggplot() +
      geom_sf(data = gulf_sf, fill = "white") +
      geom_sf(data = filter(number_sf, 
                            # pick a sequence of years to check out
                            year %in% years), 
              aes_string(col = `species`)) +
      facet_wrap(~year) +
      scale_color_viridis_c(option = "plasma") +
      labs(col = paste0("Abundance of \n",
                        str_to_lower(gsub('`','',species)))) +
             theme_void()  
  
  )
}

species_abundance_map('`WHITE HAKE`', c(1971, 1981, 1991))
