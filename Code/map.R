## Header ----
## Script name: 
##
## Purpose of script:
##
## Author: Nick Marzolf
## Date Created: 2023-03-01
## Date Modified: 
## Email: nicholas.marzolf@duke.edu
##
## load packages:  
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(spData)
library(tmap)
library(leaflet)
library(spDataLarge)
library(grid)
library(ggpubr)
##
## clear the environment if needed
rm(list = ls())
##
## set the ggplot theme
source("C:/Users/Nick Marzolf/Desktop/Research/R code/theme_nick.R")
theme_set(theme_nick())


streams_gps <- readr::read_csv('C:/Users/Nick Marzolf/Desktop/NCSU/STREAMS/La Selva GIS data/LTREB Data/GPS sites.csv')

cwd_sites <- unique(data$Stream)

cwd_coords <- streams_gps %>% 
  filter(Site %in% cwd_sites) %>% 
  st_as_sf(., coords = c('Long', 'Lat'))

# La Selva boundary
lsbs <- st_read(dsn = 'Data/Spatial/laselvaboundary.shp')

# Streams at La Selva shapefile
streams <- st_read(dsn = 'Data/Spatial/streamsclip.shp')


cr <- world %>%
  filter(name_long == 'Costa Rica')

# make a basic map of Costa Rica
cr_map <- tm_shape(cr)+   # create shape based on Costa Rica object
  tm_polygons()+          # add cr as a polygon
  tm_shape(lsbs)+         # create shape for the boundary of La Selva
  tm_dots(size = 1)       # add lsbs as a dot
cr_map


# La Selva boundary and stream layer ----

# Create a map of La Selva boundary and the stream network
map_lsbs <- tm_shape(lsbs)+            # new shape: La Selva boundary
  tm_borders()+                        # add as a border/line layer
  tm_shape(streams)+                   # new shape: stream network
  tm_lines(col = 'blue')+              # add as a line, colored blue
  tm_scale_bar(breaks = c(0, 1, 2),    # add a scale bar, with demarkations for 0, 1, and 2 km
               text.size = 0.75,       # change text size
               position = c('left',    # put the scale bar in the bottom left
                            'bottom'))
map_lsbs


map_cwd_sites <- tm_shape(lsbs)+                   # create La Selva boundary layer
  tm_borders(col = 'black')+
  tm_shape(streams)+                              # create stream network layer
  tm_lines(col = 'blue')+
  tm_shape(cwd_coords)+                             # map the locations of pH sites
  tm_symbols(size = 1,                            # change the size
             col = 'Site',                        # colored by site name
             border.col = 'black',                # with black boundary color
             palette = "viridis", n = 5)+         # change the color palette
  tm_scale_bar(breaks = c(0, 1, 2),               # add scale bar
               text.size = 0.75,
               position = c('left', 'bottom'))+
  tm_layout(inner.margins = c(.15,.01, .01, .4),  # change the margins to fit the legend and inset map
            legend.position = c('right', 'top'))+
  tm_compass(position = c('left', 'top'))         # add compass north star
map_cwd_sites

# add the inset of Costa Rica into the map
map_cwd_sites                                       # add the map with: La Selva boundary, stream layer, and points to the plot window
map_cwd_inset <- print(cr_map,                      # add the inset
                      vp = viewport(0.73, 0.205,   # change the values here to move the inset location
                                    width = 0.35,  # change the values here to change the size of the inset image
                                    height = 0.35))
