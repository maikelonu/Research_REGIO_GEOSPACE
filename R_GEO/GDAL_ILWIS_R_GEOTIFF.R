# REM Tal parece que la mejor opcion es utilizar el MOTOR interno de la version 3.8.6 de ILWIS para transformar a Geotiff
# REM dado que si preserva la georeferencia y el sistema coordenado.
# REM El problema es que genera unos archivos enormes (> 200 Mbs) y el nodata value es raro: -1e+308, tal y como se ve en QGIS
# REM La opcion es cambiarlo en R.
# REM En ILWIS entonces: export TIFF(dem_canas.mpr,dem_canas)

# Working directory is defined
setwd("/home/precis/Virtual_Folder/Export_GeoTIFF")

library(raster)
library(sp)
library(rgdal)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: barranca
# /////////////////////////////////////////////////////////////////////////////////////

dem_barranca_export <- raster("dem_barranca.tif")
Spercent_barranca_export <- raster("Spercent_barranca.tif")
uso_barranca_export <- raster("uso_barranca.tif")
clay_barranca_export <- raster("clay_barranca.tif")
density_barranca_export <- raster("density_barranca.tif")
wetness_barranca_export <- raster("wetness_barranca.tif")

NAvalue(dem_barranca_export) <- 9999
NAvalue(Spercent_barranca_export) <- 9999
NAvalue(uso_barranca_export) <- 9999
NAvalue(clay_barranca_export) <- 9999
NAvalue(density_barranca_export) <- 9999
NAvalue(wetness_barranca_export) <- 9999

dem_barranca_export[dem_barranca_export < 0] <- NA
dem_barranca_export[dem_barranca_export > 10000] <- NA

Spercent_barranca_export[Spercent_barranca_export < 0] <- NA
Spercent_barranca_export[Spercent_barranca_export > 10000] <- NA

uso_barranca_export[uso_barranca_export < 0] <- NA
uso_barranca_export[uso_barranca_export > 10000] <- NA

clay_barranca_export[clay_barranca_export < 0] <- NA
clay_barranca_export[clay_barranca_export > 10000] <- NA

density_barranca_export[density_barranca_export < 0] <- NA
density_barranca_export[density_barranca_export > 10000] <- NA

wetness_barranca_export[wetness_barranca_export < 0] <- NA
wetness_barranca_export[wetness_barranca_export > 10000] <- NA

writeRaster(dem_barranca_export, "dem_barranca_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_barranca_export, "Spercent_barranca_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_barranca_export, "uso_barranca_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_barranca_export, "clay_barranca_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_barranca_export, "density_barranca_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_barranca_export, "wetness_barranca_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: lagarto
# /////////////////////////////////////////////////////////////////////////////////////

dem_lagarto_export <- raster("dem_lagarto.tif")
Spercent_lagarto_export <- raster("Spercent_lagarto.tif")
uso_lagarto_export <- raster("uso_lagarto.tif")
clay_lagarto_export <- raster("clay_lagarto.tif")
density_lagarto_export <- raster("density_lagarto.tif")
wetness_lagarto_export <- raster("wetness_lagarto.tif")

NAvalue(dem_lagarto_export) <- 9999
NAvalue(Spercent_lagarto_export) <- 9999
NAvalue(uso_lagarto_export) <- 9999
NAvalue(clay_lagarto_export) <- 9999
NAvalue(density_lagarto_export) <- 9999
NAvalue(wetness_lagarto_export) <- 9999

dem_lagarto_export[dem_lagarto_export < 0] <- NA
dem_lagarto_export[dem_lagarto_export > 10000] <- NA

Spercent_lagarto_export[Spercent_lagarto_export < 0] <- NA
Spercent_lagarto_export[Spercent_lagarto_export > 10000] <- NA

uso_lagarto_export[uso_lagarto_export < 0] <- NA
uso_lagarto_export[uso_lagarto_export > 10000] <- NA

clay_lagarto_export[clay_lagarto_export < 0] <- NA
clay_lagarto_export[clay_lagarto_export > 10000] <- NA

density_lagarto_export[density_lagarto_export < 0] <- NA
density_lagarto_export[density_lagarto_export > 10000] <- NA

wetness_lagarto_export[wetness_lagarto_export < 0] <- NA
wetness_lagarto_export[wetness_lagarto_export > 10000] <- NA

writeRaster(dem_lagarto_export, "dem_lagarto_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_lagarto_export, "Spercent_lagarto_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_lagarto_export, "uso_lagarto_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_lagarto_export, "clay_lagarto_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_lagarto_export, "density_lagarto_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_lagarto_export, "wetness_lagarto_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: abangares
# /////////////////////////////////////////////////////////////////////////////////////

dem_abangares_export <- raster("dem_abangares.tif")
Spercent_abangares_export <- raster("Spercent_abangares.tif")
uso_abangares_export <- raster("uso_abangares.tif")
clay_abangares_export <- raster("clay_abangares.tif")
density_abangares_export <- raster("density_abangares.tif")
wetness_abangares_export <- raster("wetness_abangares.tif")

NAvalue(dem_abangares_export) <- 9999
NAvalue(Spercent_abangares_export) <- 9999
NAvalue(uso_abangares_export) <- 9999
NAvalue(clay_abangares_export) <- 9999
NAvalue(density_abangares_export) <- 9999
NAvalue(wetness_abangares_export) <- 9999

dem_abangares_export[dem_abangares_export < 0] <- NA
dem_abangares_export[dem_abangares_export > 10000] <- NA

Spercent_abangares_export[Spercent_abangares_export < 0] <- NA
Spercent_abangares_export[Spercent_abangares_export > 10000] <- NA

uso_abangares_export[uso_abangares_export < 0] <- NA
uso_abangares_export[uso_abangares_export > 10000] <- NA

clay_abangares_export[clay_abangares_export < 0] <- NA
clay_abangares_export[clay_abangares_export > 10000] <- NA

density_abangares_export[density_abangares_export < 0] <- NA
density_abangares_export[density_abangares_export > 10000] <- NA

wetness_abangares_export[wetness_abangares_export < 0] <- NA
wetness_abangares_export[wetness_abangares_export > 10000] <- NA

writeRaster(dem_abangares_export, "dem_abangares_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_abangares_export, "Spercent_abangares_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_abangares_export, "uso_abangares_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_abangares_export, "clay_abangares_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_abangares_export, "density_abangares_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_abangares_export, "wetness_abangares_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: magdalena
# /////////////////////////////////////////////////////////////////////////////////////

dem_magdalena_export <- raster("dem_magdalena.tif")
Spercent_magdalena_export <- raster("Spercent_magdalena.tif")
uso_magdalena_export <- raster("uso_magdalena.tif")
clay_magdalena_export <- raster("clay_magdalena.tif")
density_magdalena_export <- raster("density_magdalena.tif")
wetness_magdalena_export <- raster("wetness_magdalena.tif")

NAvalue(dem_magdalena_export) <- 9999
NAvalue(Spercent_magdalena_export) <- 9999
NAvalue(uso_magdalena_export) <- 9999
NAvalue(clay_magdalena_export) <- 9999
NAvalue(density_magdalena_export) <- 9999
NAvalue(wetness_magdalena_export) <- 9999

dem_magdalena_export[dem_magdalena_export < 0] <- NA
dem_magdalena_export[dem_magdalena_export > 10000] <- NA

Spercent_magdalena_export[Spercent_magdalena_export < 0] <- NA
Spercent_magdalena_export[Spercent_magdalena_export > 10000] <- NA

uso_magdalena_export[uso_magdalena_export < 0] <- NA
uso_magdalena_export[uso_magdalena_export > 10000] <- NA

clay_magdalena_export[clay_magdalena_export < 0] <- NA
clay_magdalena_export[clay_magdalena_export > 10000] <- NA

density_magdalena_export[density_magdalena_export < 0] <- NA
density_magdalena_export[density_magdalena_export > 10000] <- NA

wetness_magdalena_export[wetness_magdalena_export < 0] <- NA
wetness_magdalena_export[wetness_magdalena_export > 10000] <- NA

writeRaster(dem_magdalena_export, "dem_magdalena_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_magdalena_export, "Spercent_magdalena_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_magdalena_export, "uso_magdalena_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_magdalena_export, "clay_magdalena_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_magdalena_export, "density_magdalena_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_magdalena_export, "wetness_magdalena_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: blanco
# /////////////////////////////////////////////////////////////////////////////////////

dem_blanco_export <- raster("dem_blanco.tif")
Spercent_blanco_export <- raster("Spercent_blanco.tif")
uso_blanco_export <- raster("uso_blanco.tif")
clay_blanco_export <- raster("clay_blanco.tif")
density_blanco_export <- raster("density_blanco.tif")
wetness_blanco_export <- raster("wetness_blanco.tif")

NAvalue(dem_blanco_export) <- 9999
NAvalue(Spercent_blanco_export) <- 9999
NAvalue(uso_blanco_export) <- 9999
NAvalue(clay_blanco_export) <- 9999
NAvalue(density_blanco_export) <- 9999
NAvalue(wetness_blanco_export) <- 9999

dem_blanco_export[dem_blanco_export < 0] <- NA
dem_blanco_export[dem_blanco_export > 10000] <- NA

Spercent_blanco_export[Spercent_blanco_export < 0] <- NA
Spercent_blanco_export[Spercent_blanco_export > 10000] <- NA

uso_blanco_export[uso_blanco_export < 0] <- NA
uso_blanco_export[uso_blanco_export > 10000] <- NA

clay_blanco_export[clay_blanco_export < 0] <- NA
clay_blanco_export[clay_blanco_export > 10000] <- NA

density_blanco_export[density_blanco_export < 0] <- NA
density_blanco_export[density_blanco_export > 10000] <- NA

wetness_blanco_export[wetness_blanco_export < 0] <- NA
wetness_blanco_export[wetness_blanco_export > 10000] <- NA

writeRaster(dem_blanco_export, "dem_blanco_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_blanco_export, "Spercent_blanco_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_blanco_export, "uso_blanco_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_blanco_export, "clay_blanco_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_blanco_export, "density_blanco_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_blanco_export, "wetness_blanco_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: santarosa
# /////////////////////////////////////////////////////////////////////////////////////

dem_santarosa_export <- raster("dem_santarosa.tif")
Spercent_santarosa_export <- raster("Spercent_santarosa.tif")
uso_santarosa_export <- raster("uso_santarosa.tif")
clay_santarosa_export <- raster("clay_santarosa.tif")
density_santarosa_export <- raster("density_santarosa.tif")
wetness_santarosa_export <- raster("wetness_santarosa.tif")

NAvalue(dem_santarosa_export) <- 9999
NAvalue(Spercent_santarosa_export) <- 9999
NAvalue(uso_santarosa_export) <- 9999
NAvalue(clay_santarosa_export) <- 9999
NAvalue(density_santarosa_export) <- 9999
NAvalue(wetness_santarosa_export) <- 9999

dem_santarosa_export[dem_santarosa_export < 0] <- NA
dem_santarosa_export[dem_santarosa_export > 10000] <- NA

Spercent_santarosa_export[Spercent_santarosa_export < 0] <- NA
Spercent_santarosa_export[Spercent_santarosa_export > 10000] <- NA

uso_santarosa_export[uso_santarosa_export < 0] <- NA
uso_santarosa_export[uso_santarosa_export > 10000] <- NA

clay_santarosa_export[clay_santarosa_export < 0] <- NA
clay_santarosa_export[clay_santarosa_export > 10000] <- NA

density_santarosa_export[density_santarosa_export < 0] <- NA
density_santarosa_export[density_santarosa_export > 10000] <- NA

wetness_santarosa_export[wetness_santarosa_export < 0] <- NA
wetness_santarosa_export[wetness_santarosa_export > 10000] <- NA

writeRaster(dem_santarosa_export, "dem_santarosa_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_santarosa_export, "Spercent_santarosa_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_santarosa_export, "uso_santarosa_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_santarosa_export, "clay_santarosa_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_santarosa_export, "density_santarosa_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_santarosa_export, "wetness_santarosa_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: tenorio
# /////////////////////////////////////////////////////////////////////////////////////

dem_tenorio_export <- raster("dem_tenorio.tif")
Spercent_tenorio_export <- raster("Spercent_tenorio.tif")
uso_tenorio_export <- raster("uso_tenorio.tif")
clay_tenorio_export <- raster("clay_tenorio.tif")
density_tenorio_export <- raster("density_tenorio.tif")
wetness_tenorio_export <- raster("wetness_tenorio.tif")

NAvalue(dem_tenorio_export) <- 9999
NAvalue(Spercent_tenorio_export) <- 9999
NAvalue(uso_tenorio_export) <- 9999
NAvalue(clay_tenorio_export) <- 9999
NAvalue(density_tenorio_export) <- 9999
NAvalue(wetness_tenorio_export) <- 9999

dem_tenorio_export[dem_tenorio_export < 0] <- NA
dem_tenorio_export[dem_tenorio_export > 10000] <- NA

Spercent_tenorio_export[Spercent_tenorio_export < 0] <- NA
Spercent_tenorio_export[Spercent_tenorio_export > 10000] <- NA

uso_tenorio_export[uso_tenorio_export < 0] <- NA
uso_tenorio_export[uso_tenorio_export > 10000] <- NA

clay_tenorio_export[clay_tenorio_export < 0] <- NA
clay_tenorio_export[clay_tenorio_export > 10000] <- NA

density_tenorio_export[density_tenorio_export < 0] <- NA
density_tenorio_export[density_tenorio_export > 10000] <- NA

wetness_tenorio_export[wetness_tenorio_export < 0] <- NA
wetness_tenorio_export[wetness_tenorio_export > 10000] <- NA

writeRaster(dem_tenorio_export, "dem_tenorio_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_tenorio_export, "Spercent_tenorio_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_tenorio_export, "uso_tenorio_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_tenorio_export, "clay_tenorio_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_tenorio_export, "density_tenorio_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_tenorio_export, "wetness_tenorio_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: corobici
# /////////////////////////////////////////////////////////////////////////////////////

dem_corobici_export <- raster("dem_corobici.tif")
Spercent_corobici_export <- raster("Spercent_corobici.tif")
uso_corobici_export <- raster("uso_corobici.tif")
clay_corobici_export <- raster("clay_corobici.tif")
density_corobici_export <- raster("density_corobici.tif")
wetness_corobici_export <- raster("wetness_corobici.tif")

NAvalue(dem_corobici_export) <- 9999
NAvalue(Spercent_corobici_export) <- 9999
NAvalue(uso_corobici_export) <- 9999
NAvalue(clay_corobici_export) <- 9999
NAvalue(density_corobici_export) <- 9999
NAvalue(wetness_corobici_export) <- 9999

dem_corobici_export[dem_corobici_export < 0] <- NA
dem_corobici_export[dem_corobici_export > 10000] <- NA

Spercent_corobici_export[Spercent_corobici_export < 0] <- NA
Spercent_corobici_export[Spercent_corobici_export > 10000] <- NA

uso_corobici_export[uso_corobici_export < 0] <- NA
uso_corobici_export[uso_corobici_export > 10000] <- NA

clay_corobici_export[clay_corobici_export < 0] <- NA
clay_corobici_export[clay_corobici_export > 10000] <- NA

density_corobici_export[density_corobici_export < 0] <- NA
density_corobici_export[density_corobici_export > 10000] <- NA

wetness_corobici_export[wetness_corobici_export < 0] <- NA
wetness_corobici_export[wetness_corobici_export > 10000] <- NA

writeRaster(dem_corobici_export, "dem_corobici_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_corobici_export, "Spercent_corobici_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_corobici_export, "uso_corobici_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_corobici_export, "clay_corobici_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_corobici_export, "density_corobici_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_corobici_export, "wetness_corobici_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: salitral
# /////////////////////////////////////////////////////////////////////////////////////

dem_salitral_export <- raster("dem_salitral.tif")
Spercent_salitral_export <- raster("Spercent_salitral.tif")
uso_salitral_export <- raster("uso_salitral.tif")
clay_salitral_export <- raster("clay_salitral.tif")
density_salitral_export <- raster("density_salitral.tif")
wetness_salitral_export <- raster("wetness_salitral.tif")

NAvalue(dem_salitral_export) <- 9999
NAvalue(Spercent_salitral_export) <- 9999
NAvalue(uso_salitral_export) <- 9999
NAvalue(clay_salitral_export) <- 9999
NAvalue(density_salitral_export) <- 9999
NAvalue(wetness_salitral_export) <- 9999

dem_salitral_export[dem_salitral_export < 0] <- NA
dem_salitral_export[dem_salitral_export > 10000] <- NA

Spercent_salitral_export[Spercent_salitral_export < 0] <- NA
Spercent_salitral_export[Spercent_salitral_export > 10000] <- NA

uso_salitral_export[uso_salitral_export < 0] <- NA
uso_salitral_export[uso_salitral_export > 10000] <- NA

clay_salitral_export[clay_salitral_export < 0] <- NA
clay_salitral_export[clay_salitral_export > 10000] <- NA

density_salitral_export[density_salitral_export < 0] <- NA
density_salitral_export[density_salitral_export > 10000] <- NA

wetness_salitral_export[wetness_salitral_export < 0] <- NA
wetness_salitral_export[wetness_salitral_export > 10000] <- NA

writeRaster(dem_salitral_export, "dem_salitral_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_salitral_export, "Spercent_salitral_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_salitral_export, "uso_salitral_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_salitral_export, "clay_salitral_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_salitral_export, "density_salitral_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_salitral_export, "wetness_salitral_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: elsalto
# /////////////////////////////////////////////////////////////////////////////////////

dem_elsalto_export <- raster("dem_elsalto.tif")
Spercent_elsalto_export <- raster("Spercent_elsalto.tif")
uso_elsalto_export <- raster("uso_elsalto.tif")
clay_elsalto_export <- raster("clay_elsalto.tif")
density_elsalto_export <- raster("density_elsalto.tif")
wetness_elsalto_export <- raster("wetness_elsalto.tif")

NAvalue(dem_elsalto_export) <- 9999
NAvalue(Spercent_elsalto_export) <- 9999
NAvalue(uso_elsalto_export) <- 9999
NAvalue(clay_elsalto_export) <- 9999
NAvalue(density_elsalto_export) <- 9999
NAvalue(wetness_elsalto_export) <- 9999

dem_elsalto_export[dem_elsalto_export < 0] <- NA
dem_elsalto_export[dem_elsalto_export > 10000] <- NA

Spercent_elsalto_export[Spercent_elsalto_export < 0] <- NA
Spercent_elsalto_export[Spercent_elsalto_export > 10000] <- NA

uso_elsalto_export[uso_elsalto_export < 0] <- NA
uso_elsalto_export[uso_elsalto_export > 10000] <- NA

clay_elsalto_export[clay_elsalto_export < 0] <- NA
clay_elsalto_export[clay_elsalto_export > 10000] <- NA

density_elsalto_export[density_elsalto_export < 0] <- NA
density_elsalto_export[density_elsalto_export > 10000] <- NA

wetness_elsalto_export[wetness_elsalto_export < 0] <- NA
wetness_elsalto_export[wetness_elsalto_export > 10000] <- NA

writeRaster(dem_elsalto_export, "dem_elsalto_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_elsalto_export, "Spercent_elsalto_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_elsalto_export, "uso_elsalto_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_elsalto_export, "clay_elsalto_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_elsalto_export, "density_elsalto_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_elsalto_export, "wetness_elsalto_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: colorado
# /////////////////////////////////////////////////////////////////////////////////////

dem_colorado_export <- raster("dem_colorado.tif")
Spercent_colorado_export <- raster("Spercent_colorado.tif")
uso_colorado_export <- raster("uso_colorado.tif")
clay_colorado_export <- raster("clay_colorado.tif")
density_colorado_export <- raster("density_colorado.tif")
wetness_colorado_export <- raster("wetness_colorado.tif")

NAvalue(dem_colorado_export) <- 9999
NAvalue(Spercent_colorado_export) <- 9999
NAvalue(uso_colorado_export) <- 9999
NAvalue(clay_colorado_export) <- 9999
NAvalue(density_colorado_export) <- 9999
NAvalue(wetness_colorado_export) <- 9999

dem_colorado_export[dem_colorado_export < 0] <- NA
dem_colorado_export[dem_colorado_export > 10000] <- NA

Spercent_colorado_export[Spercent_colorado_export < 0] <- NA
Spercent_colorado_export[Spercent_colorado_export > 10000] <- NA

uso_colorado_export[uso_colorado_export < 0] <- NA
uso_colorado_export[uso_colorado_export > 10000] <- NA

clay_colorado_export[clay_colorado_export < 0] <- NA
clay_colorado_export[clay_colorado_export > 10000] <- NA

density_colorado_export[density_colorado_export < 0] <- NA
density_colorado_export[density_colorado_export > 10000] <- NA

wetness_colorado_export[wetness_colorado_export < 0] <- NA
wetness_colorado_export[wetness_colorado_export > 10000] <- NA

writeRaster(dem_colorado_export, "dem_colorado_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_colorado_export, "Spercent_colorado_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_colorado_export, "uso_colorado_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_colorado_export, "clay_colorado_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_colorado_export, "density_colorado_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_colorado_export, "wetness_colorado_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: morote
# /////////////////////////////////////////////////////////////////////////////////////

dem_morote_export <- raster("dem_morote.tif")
Spercent_morote_export <- raster("Spercent_morote.tif")
uso_morote_export <- raster("uso_morote.tif")
clay_morote_export <- raster("clay_morote.tif")
density_morote_export <- raster("density_morote.tif")
wetness_morote_export <- raster("wetness_morote.tif")

NAvalue(dem_morote_export) <- 9999
NAvalue(Spercent_morote_export) <- 9999
NAvalue(uso_morote_export) <- 9999
NAvalue(clay_morote_export) <- 9999
NAvalue(density_morote_export) <- 9999
NAvalue(wetness_morote_export) <- 9999

dem_morote_export[dem_morote_export < 0] <- NA
dem_morote_export[dem_morote_export > 10000] <- NA

Spercent_morote_export[Spercent_morote_export < 0] <- NA
Spercent_morote_export[Spercent_morote_export > 10000] <- NA

uso_morote_export[uso_morote_export < 0] <- NA
uso_morote_export[uso_morote_export > 10000] <- NA

clay_morote_export[clay_morote_export < 0] <- NA
clay_morote_export[clay_morote_export > 10000] <- NA

density_morote_export[density_morote_export < 0] <- NA
density_morote_export[density_morote_export > 10000] <- NA

wetness_morote_export[wetness_morote_export < 0] <- NA
wetness_morote_export[wetness_morote_export > 10000] <- NA

writeRaster(dem_morote_export, "dem_morote_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_morote_export, "Spercent_morote_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_morote_export, "uso_morote_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_morote_export, "clay_morote_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_morote_export, "density_morote_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_morote_export, "wetness_morote_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: canas
# /////////////////////////////////////////////////////////////////////////////////////

dem_canas_export <- raster("dem_canas.tif")
Spercent_canas_export <- raster("Spercent_canas.tif")
uso_canas_export <- raster("uso_canas.tif")
clay_canas_export <- raster("clay_canas.tif")
density_canas_export <- raster("density_canas.tif")
wetness_canas_export <- raster("wetness_canas.tif")

NAvalue(dem_canas_export) <- 9999
NAvalue(Spercent_canas_export) <- 9999
NAvalue(uso_canas_export) <- 9999
NAvalue(clay_canas_export) <- 9999
NAvalue(density_canas_export) <- 9999
NAvalue(wetness_canas_export) <- 9999

dem_canas_export[dem_canas_export < 0] <- NA
dem_canas_export[dem_canas_export > 10000] <- NA

Spercent_canas_export[Spercent_canas_export < 0] <- NA
Spercent_canas_export[Spercent_canas_export > 10000] <- NA

uso_canas_export[uso_canas_export < 0] <- NA
uso_canas_export[uso_canas_export > 10000] <- NA

clay_canas_export[clay_canas_export < 0] <- NA
clay_canas_export[clay_canas_export > 10000] <- NA

density_canas_export[density_canas_export < 0] <- NA
density_canas_export[density_canas_export > 10000] <- NA

wetness_canas_export[wetness_canas_export < 0] <- NA
wetness_canas_export[wetness_canas_export > 10000] <- NA

writeRaster(dem_canas_export, "dem_canas_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_canas_export, "Spercent_canas_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_canas_export, "uso_canas_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_canas_export, "clay_canas_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_canas_export, "density_canas_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_canas_export, "wetness_canas_export", format = "GTiff", overwrite=TRUE)

# /////////////////////////////////////////////////////////////////////////////////////
# BLOCK: tempisque
# /////////////////////////////////////////////////////////////////////////////////////

dem_tempisque_export <- raster("dem_tempisque.tif")
Spercent_tempisque_export <- raster("Spercent_tempisque.tif")
uso_tempisque_export <- raster("uso_tempisque.tif")
clay_tempisque_export <- raster("clay_tempisque.tif")
density_tempisque_export <- raster("density_tempisque.tif")
wetness_tempisque_export <- raster("wetness_tempisque.tif")

NAvalue(dem_tempisque_export) <- 9999
NAvalue(Spercent_tempisque_export) <- 9999
NAvalue(uso_tempisque_export) <- 9999
NAvalue(clay_tempisque_export) <- 9999
NAvalue(density_tempisque_export) <- 9999
NAvalue(wetness_tempisque_export) <- 9999

dem_tempisque_export[dem_tempisque_export < 0] <- NA
dem_tempisque_export[dem_tempisque_export > 10000] <- NA

Spercent_tempisque_export[Spercent_tempisque_export < 0] <- NA
Spercent_tempisque_export[Spercent_tempisque_export > 10000] <- NA

uso_tempisque_export[uso_tempisque_export < 0] <- NA
uso_tempisque_export[uso_tempisque_export > 10000] <- NA

clay_tempisque_export[clay_tempisque_export < 0] <- NA
clay_tempisque_export[clay_tempisque_export > 10000] <- NA

density_tempisque_export[density_tempisque_export < 0] <- NA
density_tempisque_export[density_tempisque_export > 10000] <- NA

wetness_tempisque_export[wetness_tempisque_export < 0] <- NA
wetness_tempisque_export[wetness_tempisque_export > 10000] <- NA

writeRaster(dem_tempisque_export, "dem_tempisque_export", format = "GTiff", overwrite=TRUE)
writeRaster(Spercent_tempisque_export, "Spercent_tempisque_export", format = "GTiff", overwrite=TRUE)
writeRaster(uso_tempisque_export, "uso_tempisque_export", format = "GTiff", overwrite=TRUE)
writeRaster(clay_tempisque_export, "clay_tempisque_export", format = "GTiff", overwrite=TRUE)
writeRaster(density_tempisque_export, "density_tempisque_export", format = "GTiff", overwrite=TRUE)
writeRaster(wetness_tempisque_export, "wetness_tempisque_export", format = "GTiff", overwrite=TRUE)


# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# END OF CODE
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////















#x <- raster("dem_canas.tif")
#x <- raster("wetness_canas.tif")

#max2 <-x@data@min

#NAvalue(x) <- 9999
#NAvalue(x) 

#x[x < 0] <- NA
#x[x > 10000] <- NA

#writeRaster(x, "pokusNA", format = "GTiff", overwrite=TRUE)
