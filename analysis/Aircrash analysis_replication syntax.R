############################################
## Spatial analysis of airplane accidents ##
############################################

# R version 3.2.3 Patched (2016-02-15 r70179)

# Clear workspace
rm(list=ls())

# Install packages
install.packages("ggplot2")
install.packages("maptools")
install.packages("raster")
install.packages("RColorBrewer")
install.packages("rgdal")
install.packages("rgeos")
install.packages("sp")
install.packages("spatstat")
install.packages("splancs")

# Load packages
library(ggplot2)
library(maptools)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(sp)
library(spatstat)
library(splancs)

# Set working directory
setwd("~/Desktop/Aircrash analysis")


############################
# Spatial data preparation #
############################

# Import aircrash data (source: NTSB aviation accident database (https://www.ntsb.gov/_layouts/ntsb.aviation/index.aspx))
events_FL <- read.table("accident_data_FL_2014.txt", sep= "|", header = T, fill = T, quote = "", row.names = NULL, stringsAsFactors = F)

# Check for duplicates
duplicated(events_FL)

# Exclude single observation with missing coordinates
events_FL <- events_FL[complete.cases(events_FL[, c(7, 8)]), ]

# Functional approach to adding missing minus signs for six observations
events_FL_corr <- events_FL 
events_FL_corr$Longitude <- with(events_FL_corr, ifelse(Longitude >= 0, -Longitude, Longitude))

# Option 1: Load shapefile (source: United States Census Bureau)
#US_shp <- readShapeSpatial("cb_2014_us_state_5m.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))

# Option 2: Download, unzip and import cartographic boundary shapefiles from USCB webpage
temp <- tempfile()
download.file("http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip", temp, mode = "w")
unzip(temp)
US_shp <- readShapeSpatial("cb_2014_us_state_5m.shp", proj4string = CRS("+proj=longlat +ellps=WGS84"))
unlink(temp)

# Crop country border to Continental US extent
US_shp_cropped <- crop(US_shp, extent(-124.848974, -66.885444, 24.396308, 49.384358)) 

# Select state border for Florida 
FL_shp <- subset(US_shp_cropped, NAME == "Florida")

# Transform SpatialPolygonsDataFrame to SpatialPolygon
FL_polygon <- SpatialPolygons(FL_shp@polygons, proj4string = FL_shp@proj4string)

# Code spatial points for subsequent analysis
FL_sp_point <- cbind(events_FL_corr$Longitude, events_FL_corr$Latitude) 
colnames(FL_sp_point) <- c("LONG", "LAT")
proj <- CRS("+proj=longlat +ellps=WGS84") 
FL_data.sp <- SpatialPointsDataFrame(coords = FL_sp_point, data = events_FL_corr, proj4string = proj)
FL_data.sp_cropped <- crop(FL_data.sp, FL_polygon) # crop to polygon
FL_df <- as.data.frame(FL_data.sp_cropped) # convert to dataframe for ggplot2


#################
# Visualization #
#################

# Plot study area
ggplot() + 
  geom_polygon(data = US_shp_cropped, aes(x = long, y = lat, group = group)) +
  geom_polygon(data = FL_shp, aes(x = long, y = lat, group = group), color = "red", fill = "red", alpha = 0.3) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  coord_equal()

# Plot study area with events
ggplot() + 
  geom_polygon(data = FL_shp, aes(x = long, y = lat, group = group)) +
  geom_point(data = FL_df, aes(x = Longitude, y = Latitude), size = 1.2, color = "red") +
  labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  coord_equal()


##########################
# Point Pattern Analysis #
##########################

## Spatial data preparation for PPA

# Change projection to equal-area projection 
FL_data.sp_aea <- (spTransform(FL_data.sp_cropped, CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))
FL_shp_aea <- (spTransform(FL_shp, CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))

# Modified function to generate ppp (source: sample code adapted from Sebastian Schutte (http://sebastianschutte.net/))
## Note: Allows for accidents occurring at the same location as they are not considered duplicates.
to_ppp <- function(events,polygon){
  events  <- as(events, "SpatialPoints")
  events  <- as(events, "ppp")
  polygon <- as(polygon, "SpatialPolygons")
  polygon <- as(polygon, "owin")
  ppp  <- ppp(events$x, events$y, window = polygon)
  return(ppp)
}

# Generate ppp object for accidents in Florida
FL_ppp <- to_ppp(FL_data.sp_aea, FL_shp_aea)

# Generate color scheme for plots
col_map <- brewer.pal(11, "RdGy")
col_map_extended <- colorRampPalette(rev(brewer.pal(11, "RdGy")))(50)


## First-order properties: Distribution of events

# Event counts 
plot(FL_ppp, pch = 21, bg = "red", cex = 0.8, alpha = 0.25, main = "Event counts")
plot(quadratcount(FL_ppp, nx = 5, ny = 5), add = T, col = "black")

# Density plot
density_FL <- density(FL_ppp, sigma = bw.diggle, edge=T)
plot(density_FL, main = "Density plot", col = col_map_extended, box = F)

# Perspective plot
persp(density(FL_ppp, sigma = bw.diggle, edge = T), ylab = "", zlab = "", xlab = "", main = "Perspective plot", ticktype = "simple", col = "red", box = T, shade = 0.75, phi = 15, theta = -70)

# Contour plot
plot(density_FL, main = "Contour of point density", col = col_map_extended, box = F)
contour(density_FL, col = "white", add = T)


## Second-order properties: Interaction between events 

n <- FL_ppp$n # extract number of events from ppp

layout(matrix(c(1,1,2,3), 2, 2, byrow = T))

# Plot observed events for visual inspection
plot(FL_shp_aea, main = "Observed points")
points(FL_ppp, pch = 21, bg = "red", cex = 0.8)

# Plot random points 
plot(FL_shp_aea, main = "Random points")
points(spsample(FL_shp_aea, n, "random"), pch = 21, bg = "red", cex = 0.8)

# Plot regular points
plot(FL_shp_aea, main = "Regular points")
points(spsample(FL_shp_aea, n, "regular"), pch = 21, bg = "red", cex = 0.8)

dev.off()

# Inhomogeneous K function (generalisation of Ripley's reduced second moment function)
plot(envelope(FL_ppp, Kinhom, correction = c("best")), main = "Inhomogeneous K function")
