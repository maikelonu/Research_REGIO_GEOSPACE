# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Construction Engineering School
# MSc.Eng. Maikel Mendez Morales
# https://www.tec.ac.cr
# Email: maikel.mendez@gmail.com; mamendez@itcr.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# Skype: maikel.mendez
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT TITLE:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INFO: This script is intended for the generation of historical IDFs curves on specific Lat/Long position
# of weather-stations for durations [5, 10, 15, 30, 60, 120, 180, 360, 720, 1440] minutes
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

# Workspace is cleared
gc(); rm(list = ls())

# Scientific notation is disabled
options(scipen=999)

# Start time is recorded
start.time <- Sys.time()

# Working directory is defined
setwd("~/Dropbox/Academics/Geospace_2025/R/climatologies")

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CRAN libraries are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
require(raster)
require(dplyr)
require(gstat)
require(maptools)
require(openxlsx)
require(readxl)
require(rgdal)
require(rgeos)
require(sp)


# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: blanco SAGA georeference is loaded
# /////////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
gc(); rm(list = ls())

# The blank map is transformed into SpatialGridDataFrame
raster.saga <- read.asciigrid("asc_blanco.asc", as.image = FALSE)

# The map is transformed into RasterLayer
raster.saga <- raster(raster.saga)

# CRTM05 projection for CR is created
CRTM05 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# CRTM05 projection is assigned to rasterbrick object
crs(raster.saga) <- CRTM05

# Object attributes are requested
# It can be seen that spa.res is different !!!!
image(raster.saga, col = topo.colors(64))

# A simple plot is requested
plot(raster.saga)

# A summery is requested
summary(raster.saga)

# A dummy raster is created
raster.saga.dummy <- raster.saga

# Dummy NA values are introduced
values(raster.saga.dummy) <- NA

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: COMPLETE 1961-1990 Costa Rica IMN Precipitation Climatology
# /////////////////////////////////////////////////////////////////////////////////////

# A saga raster brick object is loaded
clim.prec.1km <- brick("JAN_1960_DIC_1990.tif")

# A simple plot is created on raster level [10]
plot(clim.prec.1km[[10]])

# Dummy NAs values are assigned
NAvalue(clim.prec.1km) <- (-9999)

# CRTM05 projection for CR is created
CRTM05 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# CRTM05 projection is assigned to raster layer
crs(clim.prec.1km) <- CRTM05

#-----------------------------------------------------------------------------------
# Complete climatology is compiled for the period 1960-1990
#-----------------------------------------------------------------------------------

# Resampling is defined based on saga georeference
clim.prec.1km <- resample(clim.prec.1km,
                          raster.saga.dummy,
                          resample = "bilinear")

# A simple plot is created on raster level [10]
plot(clim.prec.1km[[10]])

# RasterLayers are masked based on saga Climatic-Region domain
clim.prec.1km.mask <- mask(clim.prec.1km, raster.saga)

# A simple plot is created on raster level [10]
plot(clim.prec.1km.mask[[10]])

# The mean value of the raster-brick is calculated for the entire
# period (31 years). One level raster object ONLY [mm/year]
clim.prec.1km.mean <- (calc(clim.prec.1km.mask, sum))/31

# A simple plot is created on raster level [1]
plot(clim.prec.1km.mean[[1]])

# The spatially averaged raster mean value is calculated as vector [mm/year]
clim.prec.1km.mean.vector <- as.vector(cellStats(clim.prec.1km.mean, mean, na.rm = TRUE))

# The spatially averaged raster mean value is requested [mm/year]
print(clim.prec.1km.mean.vector)

# The mean value of the raster-brick is calculated for the entire
# period (31 years) as a time-series vector object. 
# Multiple-level raster objects [mm/month]
clim.prec.1km.mean.monthly <- as.vector(cellStats(clim.prec.1km.mask, mean, na.rm = TRUE))

# The spatially averaged raster mean values are requested [mm/month]
print(clim.prec.1km.mean.monthly)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: QUASI 1961-1990 Costa Rica IMN Temperature Climatology
# /////////////////////////////////////////////////////////////////////////////////////

# A saga raster brick object is loaded
clim.temp.1km <- brick("TEMP_1960_1990_MONTHLY.tif")

# A simple plot is created on raster level [10]
plot(clim.temp.1km[[10]])

# Dummy NAs values are assigned
NAvalue(clim.temp.1km) <- (-9999)

# CRTM05 projection for CR is created
CRTM05 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# CRTM05 projection is assigned to raster layer
crs(clim.temp.1km) <- CRTM05

#-----------------------------------------------------------------------------------
# Complete climatology is compiled for the period 1960-1990
#-----------------------------------------------------------------------------------

# Resampling is defined based on saga georeference
clim.temp.1km <- resample(clim.temp.1km,
                          raster.saga.dummy,
                          resample = "bilinear")

# A simple plot is created on raster level [10]
plot(clim.temp.1km[[10]])

# RasterLayers are masked based on saga Climatic-Region domain
clim.temp.1km.mask <- mask(clim.temp.1km, raster.saga)

# A simple plot is created on raster level [10]
plot(clim.temp.1km.mask[[10]])

# The mean value of the raster-brick is calculated for the entire
# period (12 years). One level raster object ONLY [C]
clim.temp.1km.mean <- (calc(clim.temp.1km.mask, sum))/12

# A simple plot is created on raster level [1]
plot(clim.temp.1km.mean[[1]])

# The spatially averaged raster mean value is calculated as vector [C]
clim.temp.1km.mean.vector <- as.vector(cellStats(clim.temp.1km.mean, mean, na.rm = TRUE))

# The spatially averaged raster mean value is requested [C]
print(clim.temp.1km.mean.vector)

# The mean value of the raster-brick is calculated for the entire
# period (31 years) as a time-series vector object. 
# Multiple-level raster objects [C]
clim.temp.1km.mean.monthly <- as.vector(cellStats(clim.temp.1km.mask, mean, na.rm = TRUE))

# Vector is repeated 31 times to replicate entire period
clim.temp.1km.mean.monthly <- rep(clim.temp.1km.mean.monthly, 31)

# The spatially averaged raster mean values are requested [C]
print(clim.temp.1km.mean.monthly)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: QUASI 1961-1990 Costa Rica IMN ET0 Climatology
# /////////////////////////////////////////////////////////////////////////////////////

# A saga raster brick object is loaded
clim.et0.1km <- brick("ET0_1960_1990_MONTHLY.tif")

# A simple plot is created on raster level [10]
plot(clim.et0.1km[[10]])

# Dummy NAs values are assigned
NAvalue(clim.et0.1km) <- (-9999)

# CRTM05 projection for CR is created
CRTM05 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# CRTM05 projection is assigned to raster layer
crs(clim.et0.1km) <- CRTM05

#-----------------------------------------------------------------------------------
# Complete climatology is compiled for the period 1960-1990
#-----------------------------------------------------------------------------------

# Resampling is defined based on saga georeference
clim.et0.1km <- resample(clim.et0.1km,
                         raster.saga.dummy,
                         resample = "bilinear")

# A simple plot is created on raster level [10]
plot(clim.et0.1km[[10]])

# RasterLayers are masked based on saga Climatic-Region domain
clim.et0.1km.mask <- mask(clim.et0.1km, raster.saga)

# A simple plot is created on raster level [10]
plot(clim.et0.1km.mask[[10]])

# The mean value of the raster-brick is calculated for the entire
# period (12 years). One level raster object ONLY [mm/day]
clim.et0.1km.mean <- (calc(clim.et0.1km.mask, sum))/12

# A simple plot is created on raster level [1]
plot(clim.et0.1km.mean[[1]])

# The spatially averaged raster mean value is calculated as vector [mm/day]
clim.et0.1km.mean.vector <- as.vector(cellStats(clim.et0.1km.mean, mean, na.rm = TRUE))

# The spatially averaged raster mean value is requested [mm/day]
print(clim.et0.1km.mean.vector)

# The mean value of the raster-brick is calculated for the entire
# period (31 years) as a time-series vector object. 
# Multiple-level raster objects [mm/day]
clim.et0.1km.mean.monthly <- as.vector(cellStats(clim.et0.1km.mask, mean, na.rm = TRUE))

# Vector is repeated 31 times to replicate entire period
clim.et0.1km.mean.monthly <- rep(clim.et0.1km.mean.monthly, 31)

# The spatially averaged raster mean values are requested [mm/day]
print(clim.et0.1km.mean.monthly)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Entire rasterbrick object is exported to GeoTIFF format [mm/month]
writeRaster(clim.prec.1km.mask, "blanco_clim_prec_1km.tif", format="GTiff", overwrite=TRUE)

# The spatially averaged raster mean value is exported to GeoTIFF format [mm/year]
writeRaster(clim.prec.1km.mean, "blanco_clim_prec_1km_mean.tif", format="GTiff", overwrite=TRUE)

# Entire rasterbrick object is exported to GeoTIFF format [C]
writeRaster(clim.temp.1km.mask, "blanco_clim_temp_1km.tif", format="GTiff", overwrite=TRUE)

# The spatially averaged raster mean value is exported to GeoTIFF format [C]
writeRaster(clim.temp.1km.mean, "blanco_clim_temp_1km_mean.tif", format="GTiff", overwrite=TRUE)

# Entire rasterbrick object is exported to GeoTIFF format [mm/day]
writeRaster(clim.et0.1km.mask, "blanco_clim_et0_1km.tif", format="GTiff", overwrite=TRUE)

# The spatially averaged raster mean value is exported to GeoTIFF format [mm/day]
writeRaster(clim.et0.1km.mean, "blanco_clim_et0_1km_mean.tif", format="GTiff", overwrite=TRUE)

# A XLS object is created
blanco.xls.object <- list("mean_annual_prec" = clim.prec.1km.mean.vector,
                          "monthly_prec" = clim.prec.1km.mean.monthly,
                          "mean_annual_temp" = clim.temp.1km.mean.vector,
                          "monthly_temp" = clim.temp.1km.mean.monthly,
                          "mean_annual_et0" = clim.et0.1km.mean.vector,
                          "monthly_et0" = clim.et0.1km.mean.monthly)

# A XLS object is exported
write.xlsx(blanco.xls.object, file = "blanco_xls_object.xlsx", rowNames = TRUE)