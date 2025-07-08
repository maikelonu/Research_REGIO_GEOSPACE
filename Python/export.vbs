REM Tal parece que la mejor opcion es utilizar el MOTOR interno de la version 3.8.6 de ILWIS para transformar a Geotiff
REM dado que si preserva la georeferencia y el sistema coordenado.
REM El problema es que genera unos archivos enormes (> 200 Mbs) y el nodata value es raro: -1e+308, tal y como se ve en QGIS
REM La opcion es cambiarlo en R.
REM En ILWIS entonces: export TIFF(dem_canas.mpr,dem_canas)

REM barranca
REM lagarto
REM abangares
REM magdalena
REM blanco
REM santarosa
REM tenorio
REM corobici
REM salitral
REM elsalto
REM colorado
REM morote
REM canas
REM tempisque


export TIFF(dem_barranca.mpr,dem_barranca)
export TIFF(Spercent_barranca.mpr,Spercent_barranca)
export TIFF(uso_barranca.mpr,uso_barranca)
export TIFF(clay_barranca.mpr,clay_barranca)
export TIFF(density_barranca.mpr,density_barranca)
export TIFF(wetness_barranca.mpr,wetness_barranca)
export Shapefile(cat_merge_barranca.mpa,cat_merge_barranca)
export Shapefile(cat_merge_barranca.mps,cat_merge_barranca_line)

export TIFF(dem_lagarto.mpr,dem_lagarto)
export TIFF(Spercent_lagarto.mpr,Spercent_lagarto)
export TIFF(uso_lagarto.mpr,uso_lagarto)
export TIFF(clay_lagarto.mpr,clay_lagarto)
export TIFF(density_lagarto.mpr,density_lagarto)
export TIFF(wetness_lagarto.mpr,wetness_lagarto)
export Shapefile(cat_merge_lagarto.mpa,cat_merge_lagarto)
export Shapefile(cat_merge_lagarto.mps,cat_merge_lagarto_line)

export TIFF(dem_abangares.mpr,dem_abangares)
export TIFF(Spercent_abangares.mpr,Spercent_abangares)
export TIFF(uso_abangares.mpr,uso_abangares)
export TIFF(clay_abangares.mpr,clay_abangares)
export TIFF(density_abangares.mpr,density_abangares)
export TIFF(wetness_abangares.mpr,wetness_abangares)
export Shapefile(cat_merge_abangares.mpa,cat_merge_abangares)
export Shapefile(cat_merge_abangares.mps,cat_merge_abangares_line)

export TIFF(dem_magdalena.mpr,dem_magdalena)
export TIFF(Spercent_magdalena.mpr,Spercent_magdalena)
export TIFF(uso_magdalena.mpr,uso_magdalena)
export TIFF(clay_magdalena.mpr,clay_magdalena)
export TIFF(density_magdalena.mpr,density_magdalena)
export TIFF(wetness_magdalena.mpr,wetness_magdalena)
export Shapefile(cat_merge_magdalena.mpa,cat_merge_magdalena)
export Shapefile(cat_merge_magdalena.mps,cat_merge_magdalena_line)

export TIFF(dem_blanco.mpr,dem_blanco)
export TIFF(Spercent_blanco.mpr,Spercent_blanco)
export TIFF(uso_blanco.mpr,uso_blanco)
export TIFF(clay_blanco.mpr,clay_blanco)
export TIFF(density_blanco.mpr,density_blanco)
export TIFF(wetness_blanco.mpr,wetness_blanco)
export Shapefile(cat_merge_blanco.mpa,cat_merge_blanco)
export Shapefile(cat_merge_blanco.mps,cat_merge_blanco_line)

export TIFF(dem_santarosa.mpr,dem_santarosa)
export TIFF(Spercent_santarosa.mpr,Spercent_santarosa)
export TIFF(uso_santarosa.mpr,uso_santarosa)
export TIFF(clay_santarosa.mpr,clay_santarosa)
export TIFF(density_santarosa.mpr,density_santarosa)
export TIFF(wetness_santarosa.mpr,wetness_santarosa)
export Shapefile(cat_merge_santarosa.mpa,cat_merge_santarosa)
export Shapefile(cat_merge_santarosa.mps,cat_merge_santarosa_line)

export TIFF(dem_tenorio.mpr,dem_tenorio)
export TIFF(Spercent_tenorio.mpr,Spercent_tenorio)
export TIFF(uso_tenorio.mpr,uso_tenorio)
export TIFF(clay_tenorio.mpr,clay_tenorio)
export TIFF(density_tenorio.mpr,density_tenorio)
export TIFF(wetness_tenorio.mpr,wetness_tenorio)
export Shapefile(cat_merge_tenorio.mpa,cat_merge_tenorio)
export Shapefile(cat_merge_tenorio.mps,cat_merge_tenorio_line)

export TIFF(dem_carobici.mpr,dem_corobici)
export TIFF(Spercent_carobici.mpr,Spercent_corobici)
export TIFF(uso_carobici.mpr,uso_corobici)
export TIFF(clay_carobici.mpr,clay_corobici)
export TIFF(density_carobici.mpr,density_corobici)
export TIFF(wetness_carobici.mpr,wetness_corobici)
export Shapefile(cat_merge_carobici.mpa,cat_merge_corobici)
export Shapefile(cat_merge_carobici.mps,cat_merge_corobici_line)

export TIFF(dem_salitral.mpr,dem_salitral)
export TIFF(Spercent_salitral.mpr,Spercent_salitral)
export TIFF(uso_salitral.mpr,uso_salitral)
export TIFF(clay_salitral.mpr,clay_salitral)
export TIFF(density_salitral.mpr,density_salitral)
export TIFF(wetness_salitral.mpr,wetness_salitral)
export Shapefile(cat_merge_salitral.mpa,cat_merge_salitral)
export Shapefile(cat_merge_salitral.mps,cat_merge_salitral_line)

export TIFF(dem_elsalto.mpr,dem_elsalto)
export TIFF(Spercent_elsalto.mpr,Spercent_elsalto)
export TIFF(uso_elsalto.mpr,uso_elsalto)
export TIFF(clay_elsalto.mpr,clay_elsalto)
export TIFF(density_elsalto.mpr,density_elsalto)
export TIFF(wetness_elsalto.mpr,wetness_elsalto)
export Shapefile(cat_merge_elsalto.mpa,cat_merge_elsalto)
export Shapefile(cat_merge_elsalto.mps,cat_merge_elsalto_line)

export TIFF(dem_colorado.mpr,dem_colorado)
export TIFF(Spercent_colorado.mpr,Spercent_colorado)
export TIFF(uso_colorado.mpr,uso_colorado)
export TIFF(clay_colorado.mpr,clay_colorado)
export TIFF(density_colorado.mpr,density_colorado)
export TIFF(wetness_colorado.mpr,wetness_colorado)
export Shapefile(cat_merge_colorado.mpa,cat_merge_colorado)
export Shapefile(cat_merge_colorado.mps,cat_merge_colorado_line)

export TIFF(dem_morote.mpr,dem_morote)
export TIFF(Spercent_morote.mpr,Spercent_morote)
export TIFF(uso_morote.mpr,uso_morote)
export TIFF(clay_morote.mpr,clay_morote)
export TIFF(density_morote.mpr,density_morote)
export TIFF(wetness_morote.mpr,wetness_morote)
export Shapefile(cat_merge_morote.mpa,cat_merge_morote)
export Shapefile(cat_merge_morote.mps,cat_merge_morote_line)

export TIFF(dem_canas.mpr,dem_canas)
export TIFF(Spercent_canas.mpr,Spercent_canas)
export TIFF(uso_canas.mpr,uso_canas)
export TIFF(clay_canas.mpr,clay_canas)
export TIFF(density_canas.mpr,density_canas)
export TIFF(wetness_canas.mpr,wetness_canas)
export Shapefile(cat_merge_canas.mpa,cat_merge_canas)
export Shapefile(cat_merge_canas.mps,cat_merge_canas_line)

export TIFF(dem_tempisque.mpr,dem_tempisque)
export TIFF(Spercent_tempisque.mpr,Spercent_tempisque)
export TIFF(uso_tempisque.mpr,uso_tempisque)
export TIFF(clay_tempisque.mpr,clay_tempisque)
export TIFF(density_tempisque.mpr,density_tempisque)
export TIFF(wetness_tempisque.mpr,wetness_tempisque)
export Shapefile(cat_merge_tempisque.mpa,cat_merge_tempisque)
export Shapefile(cat_merge_tempisque.mps,cat_merge_tempisque_line)

export Shapefile(pc_barranca.mpp,pc_barranca)
export Shapefile(pc_lagarto.mpp,pc_lagarto)
export Shapefile(pc_abangares.mpp,pc_abangares)
export Shapefile(pc_magdalena.mpp,pc_magdalena)
export Shapefile(pc_blanco.mpp,pc_blanco)
export Shapefile(pc_santarosa.mpp,pc_santarosa)
export Shapefile(pc_tenorio.mpp,pc_tenorio)
export Shapefile(pc_salitral.mpp,pc_salitral)
export Shapefile(pc_elsalto.mpp,pc_elsalto)
export Shapefile(pc_colorado.mpp,pc_colorado)
export Shapefile(pc_morote.mpp,pc_morote)
export Shapefile(pc_canas.mpp,pc_canas)
export Shapefile(pc_tempisque.mpp,pc_tempisque)

