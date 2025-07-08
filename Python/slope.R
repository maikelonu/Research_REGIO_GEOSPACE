# Check composite spreadsheet at:
# /home/precis/Virtual_Folder/Geospace_2025/XLS/Atributos_GEO_2025.ods

# Sinuosity must be extracted from *.tbt and passed to Atributos_GEO_2025.ods
longest_magdalena.tbt

# Hydroprocessing attributes must be extracted from cat_merge_magdalena.tbt and passed to Atributos_GEO_2025.ods
# based on the concentration-point map pc_magdalena.mpp, which is "picked" directly from ILWIS GUI using the pencil icon
MapCatchmentMerge(drain_net_order.mpr,flow_dir.mpr,flow_accu.mpr,fill_base.mpr,pc_magdalena.mpp,0,1,longest_magdalena)

# This mask MUST be created prior to any other precess, 
# the domain MUST be changed to VALUE, therefore, pyramid links are BROKEN !!!
mask_magdalena = cat_merge_magdalena

# A catchment subset DEM is created
# Mean ELEV MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
dem_magdalena = iff(mask_magdalena > 0, fill_base, ?)

# A catchment X degree slope is created
fill_cut_dx_magdalena = MapFilter(dem_magdalena, dfdx)

# A catchment Y degree slope is created
fill_cut_dy_magdalena = MapFilter(dem_magdalena, dfdy)

# A catchment degree slope is created 
Sdegree_magdalena = 100 * HYP(fill_cut_dx_magdalena,fill_cut_dy_magdalena)/PIXSIZE(fill_base)

# A catchment percent slope is created
Spercent_magdalena = RADDEG(ATAN(Sdegree_magdalena/100))

# A catchment degree slope histogram is created
Spercent_magdalena.his = TableHistogram(Spercent_magdalena)

# In this case, attri_uso.mpr covers the entire domain
# A catchment landuse subset is created
# Percent forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
uso_magdalena = iff(dem_magdalena > 0, attri_uso, ?)

# In this case, soilgrids_clay.mpr covers the entire domain
# A catchment clay-content subset is created
# Clay-content MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
clay_magdalena = iff(dem_magdalena > 0,soilgrids_clay, ?)

# In this case, soilgrids_density.mpr covers the entire domain
# A catchment soil-density subset is created
# Soil-density forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
density_magdalena = iff(dem_magdalena > 0,soilgrids_density, ?)

# In this case, wetness_index.mpr covers the entire domain
# A catchment wetness-index subset is created
# Wetness-index forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
wetness_magdalena = iff(dem_magdalena > 0, wetness_index, ?)

# In this case, power_index.mpr covers the entire domain
# A catchment power-index subset is created
# power-index forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
power_magdalena = iff(dem_magdalena > 0, power_index, ?)

# In this case, sediment_index.mpr covers the entire domain
# A catchment sediment-index subset is created
# sediment-index forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
sediment_magdalena = iff(dem_magdalena > 0, sediment_index, ?)

# A catchment numeric mask is created at 30mx30m
uno_magdalena = (dem_magdalena*0) + 1

# In this case, mask_nacional.grf covers the entire Costa Rica domain at 1kmx1km
# A catchment numeric mask is created 1kmx1km
# This mask MUST be exported to ASC to be read by RStudio
asc_magdalena.mpr{dom=value.dom;vr=1.0:1.0:0.1} = MapResample(uno_magdalena,mask_nacional.grf,nearest)

# Export to ASC
export ArcInfoNAS(asc_blanco.mpr,asc_blanco)






