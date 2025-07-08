REM Check composite spreadsheet at:
REM /home/precis/Virtual_Folder/Geospace_2025/XLS/Atributos_GEO_2025.ods

REM Sinuosity must be extracted from *.tbt and passed to Atributos_GEO_2025.ods
REM longest_tempisque.tbt

REM Hydroprocessing attributes must be extracted from cat_merge_tempisque.tbt and passed to Atributos_GEO_2025.ods
REM based on the concentration-point map pc_tempisque.mpp, which is "picked" directly from ILWIS GUI using the pencil icon
REM MapCatchmentMerge(drain_net_order.mpr,flow_dir.mpr,flow_accu.mpr,fill_base_tempisque.mpr,pc_tempisque.mpp,0,1,longest_tempisque)

REM This mask MUST be created prior to any other precess, 
REM the domain MUST be changed to VALUE, therefore, pyramid links are BROKEN !!!
REM mask_tempisque = cat_merge_tempisque

REM A catchment subset DEM is created
REM Mean ELEV MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
dem_tempisque = iff(mask_tempisque > 0, fill_base_tempisque, ?)

REM A catchment X degree slope is created
fill_cut_dx_tempisque = MapFilter(dem_tempisque, dfdx)

REM A catchment Y degree slope is created
fill_cut_dy_tempisque = MapFilter(dem_tempisque, dfdy)

REM A catchment degree slope is created 
Sdegree_tempisque = 100 * HYP(fill_cut_dx_tempisque,fill_cut_dy_tempisque)/PIXSIZE(fill_base_tempisque)

REM A catchment percent slope is created
Spercent_tempisque = RADDEG(ATAN(Sdegree_tempisque/100))

REM A catchment degree slope histogram is created
REM Spercent_tempisque.his = TableHistogram(Spercent_tempisque)

REM In this case, attri_uso.mpr covers the entire domain
REM A catchment landuse subset is created
REM Percent forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
uso_tempisque = iff(dem_tempisque > 0, uso_1992_gte, ?)

REM In this case, soilgrids_clay.mpr covers the entire domain
REM A catchment clay-content subset is created
REM Clay-content MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
clay_tempisque = iff(dem_tempisque > 0,soil_clay_gte, ?)

REM In this case, soilgrids_density.mpr covers the entire domain
REM A catchment soil-density subset is created
REM Soil-density forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
density_tempisque = iff(dem_tempisque > 0,soil_density_gte, ?)

REM In this case, wetness_index.mpr covers the entire domain
REM A catchment wetness-index subset is created
REM Wetness-index forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
wetness_tempisque = iff(dem_tempisque > 0, wetness_index_tempisque, ?)

REM In this case, power_index.mpr covers the entire domain
REM A catchment power-index subset is created
REM power-index forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
REM power_tempisque = iff(dem_tempisque > 0, power_index, ?)

REM In this case, sediment_index.mpr covers the entire domain
REM A catchment sediment-index subset is created
REM sediment-index forest MUST be extracted from histogram and passed to Atributos_GEO_2025.ods
REM sediment_tempisque = iff(dem_tempisque > 0, sediment_index, ?)

REM A catchment numeric mask is created at 30mx30m
uno_tempisque = (dem_tempisque*0) + 1

REM In this case, mask_nacional.grf covers the entire Costa Rica domain at 1kmx1km
REM A catchment numeric mask is created 1kmx1km
REM This mask MUST be exported to ASC to be read by RStudio
asc_tempisque.mpr{dom=value.dom;vr=1.0:1.0:0.1} = MapResample(uno_tempisque,mask_nacional.grf,nearest)

REM Export to ASC
export ArcInfoNAS(asc_tempisque.mpr,asc_tempisque)



