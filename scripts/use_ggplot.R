# Script to explore the data using sf and ggplot

library(sf)
library(tidyverse)

# load data
explan <- readRDS("data/explan.RDS")
number <- readRDS("data/number.RDS")
gulf <- readRDS("data/gulf.RDS")

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
ggplot() +
  geom_sf(data = gulf_sf, fill = "white") +
  geom_sf(data = filter(explan_sf, 
                        # pick a sequence of years to check out
                        year %in% c(1971, 1981, 1991, 2001, 2011, 2019)),
          aes(col = bottom.temperature)) +
  facet_wrap(~year) +
  scale_color_viridis_c(option = "viridis") +
  labs(color = "Bottom \ntemperature") +
  theme_void() 


## ggplot map of abundance for a given species

ggplot() +
  geom_sf(data = gulf_sf, fill = "white") +
  geom_sf(data = filter(number_sf, 
                        # pick a sequence of years to check out
                        year %in% c(1971, 1981, 1991, 2001, 2011, 2019)), 
          aes(col = `WHITE HAKE`)) +
  facet_wrap(~year) +
  scale_color_viridis_c(option = "plasma") +
  labs(col = "Abundance of \nWhite hake") +
  theme_void() 
