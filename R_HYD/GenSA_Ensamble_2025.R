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
# /////////////////////////////////////////////////////////// /////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT TITLE:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# To be defined
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INFO: To be defined
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

# Start the clock!
ptm <- proc.time()

# Working directory is defined
setwd("~/Dropbox/Academics/IDF_CC_tool_CANADA/R_scripts/GR2M_REGIO")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: CRAN libraries are loaded
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
require(airGR)
require(airGRteaching)
require(DEoptim)
require(DescTools)
require(dplyr)
require(DWBmodelUN)
require(FME)
require(forecast)
require(GenSA)
require(ggplot2)
require(hydroGOF)
require(hydromad)
require(hydrosystems)
require(lubridate)
require(matrixStats)
require(outliers)
require(pastecs)
require(plyr)
require(reshape)
require(reshape2)
require(tidyr)
require(viridis)
require(openxlsx)
require(weathermetrics)
require(xts)
require(zoo)
require(scales)
require(clusterSim)
require(ggradar)
require(fmsb)

# In case these libraries need to be installed:
# devtools::install_github('tanerumit/hydrosystems')
# devtools::install_github("JosephGuillaume/hydromad")
# devtools::install_github("dazamora/DWBmodelUN")
# devtools::install_github("ricardo-bion/ggradar")
# install.packages("~/Downloads/rgdal_1.6-7.tar.gz", repos=NULL, type="source")
# install.packages("~/Downloads/GenSA_1.1.15.tar.gz", repos=NULL, type="source")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Weather Station Metadata + Input Data
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Weather station ETCCDI indices are loaded
df.ens.metrics <- read.delim("ENS_COMPLEX02.txt", header = TRUE, sep = "\t")

df.ens.metrics <- subset(df.ens.metrics, model != "ENSEMBLE")

df.ens.metrics <- subset(df.ens.metrics, model != "ABCD")

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: GenSA
#--------------------------------------------------------------------------------------------------------------

df.ens.metrics.gensa <- df.ens.metrics[,c(1,10,12)]
df.ens.metrics.gensa <- subset(df.ens.metrics.gensa, stage == "CAL")

df.ens.metrics.gensa.gr2m <- subset(df.ens.metrics.gensa, model == "GR2M")
min(df.ens.metrics.gensa.gr2m$kge)
max(df.ens.metrics.gensa.gr2m$kge)
mean(df.ens.metrics.gensa.gr2m$kge)
sd(df.ens.metrics.gensa.gr2m$kge)

df.ens.metrics.gensa.budyko <- subset(df.ens.metrics.gensa, model == "BUDYKO")
min(df.ens.metrics.gensa.budyko$kge)
max(df.ens.metrics.gensa.budyko$kge)
mean(df.ens.metrics.gensa.budyko$kge)
sd(df.ens.metrics.gensa.budyko$kge)

df.ens.metrics.gensa.awbd <- subset(df.ens.metrics.gensa, model == "AWBM")
min(df.ens.metrics.gensa.awbd$kge)
max(df.ens.metrics.gensa.awbd$kge)
mean(df.ens.metrics.gensa.awbd$kge)
sd(df.ens.metrics.gensa.awbd$kge)

df.ens.metrics.gensa.ihacres <- subset(df.ens.metrics.gensa, model == "IHACRES")
min(df.ens.metrics.gensa.ihacres$kge)
max(df.ens.metrics.gensa.ihacres$kge)
mean(df.ens.metrics.gensa.ihacres$kge)
sd(df.ens.metrics.gensa.ihacres$kge)

df.ens.metrics.gensa.simhyd <- subset(df.ens.metrics.gensa, model == "SIMHYD")
min(df.ens.metrics.gensa.simhyd$kge)
max(df.ens.metrics.gensa.simhyd$kge)
mean(df.ens.metrics.gensa.simhyd$kge)
sd(df.ens.metrics.gensa.simhyd$kge)
# 
# df.ens.metrics.gensa.abcd <- subset(df.ens.metrics.gensa, model == "abcd")
# min(df.ens.metrics.gensa.abcd$kge)
# max(df.ens.metrics.gensa.abcd$kge)
# mean(df.ens.metrics.gensa.abcd$kge)
# sd(df.ens.metrics.gensa.abcd$kge)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Normalization
#--------------------------------------------------------------------------------------------------------------

df.ens.metrics.normalization <- df.ens.metrics

# Test normalization
df.ens.metrics.normalization$kge.norm <- data.Normalization(df.ens.metrics.normalization$kge, type="n6", normalization="column")
df.ens.metrics.normalization$nse.norm <- data.Normalization(df.ens.metrics.normalization$nse, type="n4", normalization="column")
df.ens.metrics.normalization$r2.norm <- data.Normalization(df.ens.metrics.normalization$corr, type="n4", normalization="column")
df.ens.metrics.normalization$d.norm <- data.Normalization(df.ens.metrics.normalization$d, type="n4", normalization="column")
df.ens.metrics.normalization$md.norm <- data.Normalization(df.ens.metrics.normalization$md, type="n4", normalization="column")

# rsr
df.ens.metrics.normalization$rsr.norm <- 1 - df.ens.metrics.normalization$rsr
df.ens.metrics.normalization$rsr.norm <- data.Normalization(df.ens.metrics.normalization$rsr.norm, type="n4", normalization="column")

# pbias
df.ens.metrics.normalization$pbias.norm <- abs((df.ens.metrics.normalization$pbias)/100)
df.ens.metrics.normalization$pbias.norm <- 1- df.ens.metrics.normalization$pbias.norm
df.ens.metrics.normalization$pbias.norm <- data.Normalization(df.ens.metrics.normalization$pbias.norm, type="n4", normalization="column")

df.ens.metrics.normalization$uni <- (df.ens.metrics.normalization$kge.norm +
                                       #df.ens.metrics.normalization$nse.norm +
                                       #df.ens.metrics.normalization$r2.norm +
                                       df.ens.metrics.normalization$d.norm  + 
                                       df.ens.metrics.normalization$rsr.norm + 
                                       df.ens.metrics.normalization$pbias.norm)

#df.ens.metrics.normalization$uni <- rescale(df.ens.metrics.normalization$uni, to=c(0,1))

df.ens.metrics.normalization$ranks <- rank(df.ens.metrics.normalization$uni)

df.ens.metrics.normalization$ranks <- 5 - df.ens.metrics.normalization$ranks

df.ens.metrics.normalization$uni02 <- (df.ens.metrics.normalization$kge +
                                         #df.ens.metrics.normalization$nse +
                                         #df.ens.metrics.normalization$r2.norm +
                                         df.ens.metrics.normalization$d  + 
                                         (1 - df.ens.metrics.normalization$rsr) + 
                                         (1-abs((df.ens.metrics.normalization$pbias)/100)))

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Radar Charts
#--------------------------------------------------------------------------------------------------------------

# A new scale color is defined
scale01 <- c("black", "blue")
scale01 <- c("#a6761d","#ca0020","#4dac26", "#2166ac", "#f4a582", "#0571b0")  
scale01 <- c("#377eb8","#e41a1c", "#4daf4a", "black", "#ff7f00", "blue", "orange")
scale01 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
scale01 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#666666")
scale02 <- c("#a6761d","#ca0020","#4dac26", "#2166ac", "#0571b0")  

#-----------------------------------------------------------------
# Heatmap KGE
#-----------------------------------------------------------------

# Radar variables are selected
df.ens.radar <- df.ens.metrics.normalization[,c(10,12,13,1)]

# Reorder factor levels for 'category'
df.ens.radar$catchment <- factor(df.ens.radar$catchment, levels = c("TNR", "TMP", "SNT", "SLR", "MRT", "LAG", "EST", "CRB", "CNS", "BLC", "BAR", "ABG"))

# Radar variables are rounded
df.ens.radar$uni02 <- round(df.ens.radar$kge,2)

ggplot() +
  geom_tile(aes(x = model,y = catchment,fill = kge),data=df.ens.radar,size = 0.25) +
  facet_grid(rows = stage ~ .) +
  geom_text(aes(x = model,y = catchment,label = uni02),data=df.ens.radar,
            colour = '#333333',size = 4.0,parse = FALSE) +
  scale_fill_gradient2(
    low = "#ca0020",    # Color for low values
    mid = "#ffffb3",       # Color for the midpoint
    high = "#1f78b4",   # Data value corresponding to the 'mid' color
    midpoint =  0.70) +
  xlab("Model") +
  ylab("Catchment") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#-----------------------------------------------------------------
# Heatmap PBIAS
#-----------------------------------------------------------------

# Radar variables are selected
df.ens.radar <- df.ens.metrics.normalization[,c(10,12,13,3)]

# Reorder factor levels for 'category'
df.ens.radar$catchment <- factor(df.ens.radar$catchment, levels = c("TNR", "TMP", "SNT", "SLR", "MRT", "LAG", "EST", "CRB", "CNS", "BLC", "BAR", "ABG"))

# Radar variables are rounded
df.ens.radar$uni02 <- round(df.ens.radar$pbias,1)

ggplot() +
  geom_tile(aes(x = model,y = catchment,fill = pbias),data=df.ens.radar,size = 0.25) +
  facet_grid(rows = stage ~ .) +
  geom_text(aes(x = model,y = catchment,label = uni02),data=df.ens.radar,
            colour = '#333333',size = 4.0,parse = FALSE) +
  scale_fill_gradient2(
    low = "#ca0020",    # Color for low values
    mid = "#ffffb3",       # Color for the midpoint
    high = "#1f78b4",   # Data value corresponding to the 'mid' color
    midpoint =  0.0) +
  xlab("Model") +
  ylab("Catchment") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#-----------------------------------------------------------------
# Heatmap D
#-----------------------------------------------------------------

# Radar variables are selected
df.ens.radar <- df.ens.metrics.normalization[,c(10,12,13,5)]

# Reorder factor levels for 'category'
df.ens.radar$catchment <- factor(df.ens.radar$catchment, levels = c("TNR", "TMP", "SNT", "SLR", "MRT", "LAG", "EST", "CRB", "CNS", "BLC", "BAR", "ABG"))

# Radar variables are rounded
df.ens.radar$uni02 <- round(df.ens.radar$d,2)

ggplot() +
  geom_tile(aes(x = model,y = catchment,fill = d),data=df.ens.radar,size = 0.25) +
  facet_grid(rows = stage ~ .) +
  geom_text(aes(x = model,y = catchment,label = uni02),data=df.ens.radar,
            colour = '#333333',size = 4.0,parse = FALSE) +
  scale_fill_gradient2(
    low = "#ca0020",    # Color for low values
    mid = "#ffffb3",       # Color for the midpoint
    high = "#1f78b4",   # Data value corresponding to the 'mid' color
    midpoint =  0.85) +
  xlab("Model") +
  ylab("Catchment") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#-----------------------------------------------------------------
# Heatmap RSR
#-----------------------------------------------------------------

# Radar variables are selected
df.ens.radar <- df.ens.metrics.normalization[,c(10,12,13,9)]

# Reorder factor levels for 'category'
df.ens.radar$catchment <- factor(df.ens.radar$catchment, levels = c("TNR", "TMP", "SNT", "SLR", "MRT", "LAG", "EST", "CRB", "CNS", "BLC", "BAR", "ABG"))

# Radar variables are rounded
df.ens.radar$uni02 <- round(df.ens.radar$rsr,2)

ggplot() +
  geom_tile(aes(x = model,y = catchment,fill = rsr),data=df.ens.radar,size = 0.25) +
  facet_grid(rows = stage ~ .) +
  geom_text(aes(x = model,y = catchment,label = uni02),data=df.ens.radar,
            colour = '#333333',size = 4.0,parse = FALSE) +
  scale_fill_gradient2(
    low = "#1f78b4",    # Color for low values
    mid = "#ffffb3",       # Color for the midpoint
    high = "#ca0020",   # Data value corresponding to the 'mid' color
    midpoint =  0.65) +
  xlab("Model") +
  ylab("Catchment") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# An export/import data.frame is created to free RA memory
save(df.ens.radar, file = "windows.RData")

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Radar Calibration
#--------------------------------------------------------------------------------------------------------------

df.ens.radar.cal <- subset(df.ens.radar, stage == "CAL")

df.wide.cal <- dcast(data = df.ens.radar.cal,
                     formula = model ~ catchment,
                     fun.aggregate = mean,
                     value.var = "uni02")

ggradar::ggradar(df.wide.cal)

# v.temp_kappa_HadGEM2_RegCM47
ggradar::ggradar(df.wide.cal,
                 values.radar = c(0.0,0.50,1),
                 #grid.min = 0.50,
                 #grid.mid = 0.75,
                 #grid.max = 1.0,
                 background.circle.colour = "white",
                 axis.line.colour = "#636363",
                 gridline.min.colour = "#636363",
                 gridline.mid.colour = "#636363",
                 gridline.max.colour = "#636363",
                 label.gridline.max = FALSE,
                 grid.line.width =2.5,
                 group.line.width = 1.25,
                 group.point.size = 3,
                 grid.label.size = 10,
                 background.circle.transparency = 0.0,
                 gridline.min.linetype = 2,
                 gridline.mid.linetype = 2,
                 gridline.max.linetype = 1,
                 legend.title = "UNI-Composite-Score",
                 legend.position = "right",
                 legend.text.size = 18,
                 axis.label.size = 8, 
                 group.colours = scale02)

# --------------------------------------------

rownames(df.wide.cal) <- c("AWBM", "BUDYKO", "GR2M", "IHACRES", "SIMHYD" )

df.gg.cal <- df.wide.cal[,-c(1)]

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
df.gg.cal <- rbind(rep(3.63,12) , rep(2.93,12) , df.gg.cal)

# plot with default options:
radarchartcirc(df.gg.cal, axistype=1 , pcol=scale02, maxmin=TRUE,
               #custom polygon
               #pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
               #custom the grid
               cglcol="#636363", axislabcol="#636363",
               caxislabels=seq(2.93,3.63,0.20),
               cglwd=0.8,
               plwd=2.5,
               cglty=1,
               plty=c(1,1,1,1,1),
               vlcex=0.8 )

# Add a legend
legend(x=1.5, y=1, legend = rownames(df.gg.cal[-c(1,2),]),
       bty = "n", pch=20 , col=scale02 , text.col = "grey", cex=1.2, pt.cex=3)
# --------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Radar validation
#--------------------------------------------------------------------------------------------------------------

df.ens.radar.val <- subset(df.ens.radar, stage == "VAL")

df.wide.val <- dcast(data = df.ens.radar.val,
                     formula = model ~ catchment,
                     fun.aggregate = mean,
                     value.var = "uni02")

ggradar::ggradar(df.wide.val)

# v.temp_kappa_HadGEM2_RegCM47
ggradar::ggradar(df.wide.val,
                 #values.radar = c(0.4,0.5,1),
                 grid.min = 0,
                 grid.mid = 3,
                 grid.max = 4.6,
                 background.circle.colour = "white",
                 axis.line.colour = "#636363",
                 gridline.min.colour = "#636363",
                 gridline.mid.colour = "#636363",
                 gridline.max.colour = "#636363",
                 label.gridline.max = FALSE,
                 grid.line.width =2.5,
                 group.line.width = 1.25,
                 group.point.size = 3,
                 grid.label.size = 10,
                 background.circle.transparency = 0.0,
                 gridline.min.linetype = 2,
                 gridline.mid.linetype = 2,
                 gridline.max.linetype = 1,
                 legend.title = "UNI-Composite-Score",
                 legend.position = "right",
                 legend.text.size = 18,
                 axis.label.size = 8, 
                 group.colours = scale02)

# --------------------------------------------

rownames(df.wide.val) <- c("AWBM", "BUDYKO", "GR2M", "IHACRES", "SIMHYD" )

df.gg.val <- df.wide.val[,-c(1)]

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
df.gg.val <- rbind(rep(3.44,12) , rep(1.97,12) , df.gg.val)

# plot with default options:
radarchartcirc(df.gg.val, axistype=1 , pcol=scale02, maxmin=TRUE,
               #custom polygon
               #pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
               #custom the grid
               cglcol="#636363", axislabcol="#636363",
               caxislabels=seq(1.97,3.44,0.20),
               cglwd=0.8,
               plwd=2.5,
               cglty=1,
               plty=c(1,1,1,1,1),
               vlcex=0.8 )

# Add a legend
legend(x=1.5, y=1, legend = rownames(df.gg.val[-c(1,2),]),
       bty = "n", pch=20 , col=scale02 , text.col = "grey", cex=1.2, pt.cex=3)
# --------------------------------------------


#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: General Boxplots
#--------------------------------------------------------------------------------------------------------------

df.aggregate <- aggregate(x=df.ens.metrics$kge,
                          by=list(df.ens.metrics$model,
                                  df.ens.metrics$par_number,
                                  df.ens.metrics$comp_number,
                                  df.ens.metrics$stage),
                          FUN=mean)

names(df.aggregate) <- c("model", "par_number", "comp_number","stage","kge")



# An export/import data.frame is created to free RA memory
save(df.ens.metrics, df.aggregate, file = "windows.RData")

ggplot() +
  geom_boxplot(aes(y = kge,
                   x = as.factor(par_number),
                   colour = stage,
                   fill = model,
                   # CRITICAL: Define a single, explicit grouping variable for dodging
                   group = interaction(par_number, stage, model)),
               data = df.ens.metrics,
               outlier.shape = 3,
               # Ensure the width is 0.75
               position = position_dodge(width = 0.75)) +
  geom_point(aes(x = as.factor(par_number),
                 y = kge,
                 colour = stage,
                 # CRITICAL: Use the same explicit grouping variable for the points
                 group = interaction(par_number, stage, model)),
             data = df.aggregate,
             # Ensure the width is 0.75
             position = position_dodge(width = 0.75)) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  geom_hline(data=df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  ggtitle("KGE") +
  xlab("Number of Parameters") +
  ylab("KGE") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#KGE Performance with number of parameters
ggplot() +
  geom_boxplot(aes(y = kge,
                   x = par_number,
                   colour = stage, 
                   fill = as.factor(model)),
               data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12,
               position = position_dodge(preserve = "single")) +
    scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7.0,min.n = 7.0)) +
  stat_summary(aes(x = par_number,y = kge,
                   shape = stage,
                   colour = model,
                   group = stage,
                   size = 1.75,),
               data=df.ens.metrics,size = 4.0,
               fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',
               position = position_dodge(width = 0.50)) +
    geom_hline(data=df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  ggtitle("KGE") +
  xlab("Number of Parameters") +
  ylab("KGE") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#KGE Performance with number of storages
ggplot() +
  geom_boxplot(aes(y = kge,
                   x = comp_number,
                   colour = stage, 
                   fill = as.factor(model)),
               data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5.0,min.n = 5.0)) +
  stat_summary(aes(x = comp_number,y = kge,
                   shape = stage,
                   colour = model,
                   size = 1.75,),
               data=df.ens.metrics,size = 4.0,
               fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75)) +
  
  geom_hline(data=df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  ggtitle("KGE") +
  xlab("Number of Storages") +
  ylab("KGE") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#KGE
ggplot() +
  geom_boxplot(aes(y = kge,x = model,colour = stage, fill = as.factor(par_number)),data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = kge,shape = stage,colour = stage),data=df.ens.metrics,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75)) +
  geom_hline(data=df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  ggtitle("KGE") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#NSE
ggplot() +
  geom_boxplot(aes(y = nse,x = model,colour = stage, fill = as.factor(par_number)),data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = nse,shape = stage,colour = stage),data=df.ens.metrics,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75))+
  geom_hline(data=df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2)+
  ggtitle("NSE") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#PBIAS
ggplot() +
  geom_boxplot(aes(y = pbias,x = model,colour = stage, fill = as.factor(par_number)),data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = pbias,shape = stage,colour = stage),data=df.ens.metrics,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75)) +
  geom_hline(data=df.ens.metrics,yintercept = 15,size = 0.555,linetype = 2) +
  geom_hline(data=df.ens.metrics,yintercept = -15,size = 0.555,linetype = 2) +
  geom_hline(data=df.ens.metrics,yintercept = 0.0,size = 0.555,linetype = 2) +
  ggtitle("PBIAS") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#R2
ggplot() +
  geom_boxplot(aes(y = corr,x = model,colour = stage, fill = as.factor(par_number) ),data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = corr,shape = stage,colour = stage),data=df.ens.metrics,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75))+
  geom_hline(data=df.ens.metrics,yintercept = 0.60,size = 0.555,linetype = 2)+
  ggtitle("R2") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#D
ggplot() +
  geom_boxplot(aes(y = d,x = model,colour = stage, fill = as.factor(par_number)),data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = d,shape = stage,colour = stage),data=df.ens.metrics,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75))+
  geom_hline(data=df.ens.metrics,yintercept = 0.75,size = 0.555,linetype = 2)+
  ggtitle("D") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#rsr
ggplot() +
  geom_boxplot(aes(y = rsr,x = model,colour = stage, fill = as.factor(par_number)),data=df.ens.metrics,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = rsr,shape = stage,colour = stage),data=df.ens.metrics,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75))+
  geom_hline(data=df.ens.metrics,yintercept = 0.70,size = 0.555,linetype = 2)+
  ggtitle("RSR") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#UNI2
ggplot() +
  geom_boxplot(aes(y = uni02,x = model,colour = stage, fill = as.factor(par_number)),data=df.ens.metrics.normalization,size = 0.75, outlier.size = 4, outlier.shape = 12) +
  scale_color_manual(values = scale01) +
  theme_bw() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  stat_summary(aes(x = model,y = uni02,shape = stage,colour = stage),df.ens.metrics.normalization,size = 4.0,fun.data = mean_sdl,fun.args = list(mult = 1),geom = 'point',position = position_dodge(width = 0.75)) +
  #geom_hline(data=df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  ggtitle("UNI-Composite-Score") +
  xlab("MODEL") +
  ylab("OF") +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Mean Error-Plots
#--------------------------------------------------------------------------------------------------------------

# KGE performance vs parameters
ggplot() +
  stat_summary(aes(x = par_number, y = kge,
                   colour = stage,
                   shape = model),
               data = df.ens.metrics,
               size = 1.5,
               linewidth = 1.0,
               fun.data = mean_cl_boot,
               # Key change: Add position_dodge()
               position = position_dodge(width = 0.50)) +
  geom_hline(data= df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  xlab("Number of Parameters") +
  ylab("KGE") +
  scale_shape_manual(values = c(15, 16, 17, 18, 8)) +
  scale_color_manual(values = c("#2166ac", "#ca0020")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7.0,min.n = 7.0)) +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif")) +
  theme_bw()

# KGE performance vs storages
ggplot() +
  stat_summary(aes(x = comp_number, y = kge,
                   colour = stage,
                   shape = model),
               data = df.ens.metrics,
               size = 1.5,
               linewidth = 1.0,
               fun.data = mean_cl_boot,
               # Key change: Add position_dodge()
               position = position_dodge(width = 0.50)) +
  geom_hline(data= df.ens.metrics,yintercept = 0.50,size = 0.555,linetype = 2) +
  xlab("Number of Storages") +
  ylab("KGE") +
  scale_shape_manual(values = c(15, 16, 17, 18, 8)) +
  scale_color_manual(values = c("#2166ac", "#ca0020")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7.0,min.n = 7.0)) +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif")) +
  theme_bw()


ggplot() +
  stat_summary(aes(x = par_number,y = nse,colour = stage,shape = model),data=df.ens.metrics,size = 0.99,fun.data = mean_sdl,fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw()

ggplot() +
  stat_summary(aes(x = par_number,y = pbias,colour = stage,shape = model),data=df.ens.metrics,size = 0.99,fun.data = mean_sdl,fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw()

ggplot() +
  stat_summary(aes(x = par_number,y = d,colour = stage,shape = model),data=df.ens.metrics,size = 0.99,fun.data = mean_sdl,fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw()

ggplot() +
  stat_summary(aes(x = par_number,y = rsr,colour = stage,shape = model),data=df.ens.metrics,size = 0.99,fun.data = mean_sdl,fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw()

ggplot() +
  stat_summary(aes(x = par_number,y = uni02,colour = stage,shape = model),data=df.ens.metrics.normalization,size = 0.99,fun.data = mean_sdl,fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  theme_bw()

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: KGE-UNI02 performance vs duration period
#--------------------------------------------------------------------------------------------------------------

# Duration of calibration and validation periods are added
df.ens.metrics.normalization$dur_cal <- 0

counter01 <- length(df.ens.metrics.normalization$kge)

for(i in 1:counter01) {
  
  if(df.ens.metrics.normalization$catchment[i] == "TMP" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 5.5 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "SNT" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 5.75 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "SLR" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 7.08 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "MRT" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 8.33 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "LAG" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 11.67 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "EST" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 7.08 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "CRB" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 9.08 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "CNS" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 14.08 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "BLC" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 6.75 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "BAR" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 6.08 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "ABG" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 13.33 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "TNR" & df.ens.metrics.normalization$stage[i] == "CAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 11.25 }

}

# A disaggregation baseflow loop is executed to separate baseflow by month
for(i in 1:counter01) {
  
  if(df.ens.metrics.normalization$catchment[i] == "TMP" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 3.17 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "SNT" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 3.08 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "SLR" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 3.67 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "MRT" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 4.50 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "LAG" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 6.42 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "EST" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 4.00 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "CRB" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 9.92 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "CNS" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 9.00 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "BLC" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 4.00 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "BAR" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 7.58 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "ABG" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 6.42 }
  
  else if(df.ens.metrics.normalization$catchment[i] == "TNR" & df.ens.metrics.normalization$stage[i] == "VAL") {
    df.ens.metrics.normalization$dur_cal[i] <- 5.50 }
  
}

ggplot() +
  geom_point(aes(x = dur_cal, y = kge, colour = stage, shape = catchment), data = df.ens.metrics.normalization, size = 3.99) +
  geom_smooth(aes(x = dur_cal, y = kge, colour = stage, fill = stage), data = df.ens.metrics.normalization, method = lm, alpha = 0.40) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0, min.n = 10.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0, min.n = 10.0)) +
  facet_wrap(facets = ~ model) +
  scale_shape_manual(values = c(12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)) +
  # --- Manual Color Selection ---
  scale_colour_manual(values = c("VAL" = "#ca0020", "CAL" = "#2166ac")) +
  scale_fill_manual(values = c("VAL" = "#ca0020", "CAL" = "#2166ac")) +
  xlab("Duration Period (Years)") +
  ylab("KGE") +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), text=element_text(size=14, family="serif"))

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Correlation Q vs. Prec.
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Weather station ETCCDI indices are loaded
df.ens.corr <- read.delim("ENS_CORR.txt", header = TRUE, sep = "\t")

# observed
cor.test(df.ens.corr$observed_Q, df.ens.corr$observed_P)

# gr2m
cor.test(df.ens.corr$Q_gr2m, df.ens.corr$observed_P)

# abcd
cor.test(df.ens.corr$Q_abcd, df.ens.corr$observed_P)

# budyko
cor.test(df.ens.corr$Q.budyko, df.ens.corr$observed_P)

# awbm
cor.test(df.ens.corr$Q_awbm, df.ens.corr$observed_P)

# ihacres
cor.test(df.ens.corr$Q_ihacres, df.ens.corr$observed_P)

# simhyd
cor.test(df.ens.corr$Q.simhyd, df.ens.corr$observed_P)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: FDC Flow Duration Curves and Hydrographs
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: canas
#--------------------------------------------------------------------------------------------------------------

#----------- canas ----------
workbook <- loadWorkbook("Multi_Hydro_Ensemble_canas_SECOND-ROUND.xlsx")
sheet_name.cal <- "df.ensemble.out.cal" # Replace with the name of the sheet you want to export
data_from_sheet.cal <- read.xlsx(workbook, sheet = sheet_name.cal)
sheet_name.val <- "df.ensemble.out.val" # Replace with the name of the sheet you want to export
data_from_sheet.val <- read.xlsx(workbook, sheet = sheet_name.val)
data_from_sheet.cal$catchment <- c("canas")
data_from_sheet.val$catchment <- c("canas")
data_from_sheet.cal$stage <- c("CAL")
data_from_sheet.val$stage <- c("VAL")
data_from_sheet.val[c(1:32), 13] <- c("CAL")

# abcd and mean models are removed
data_from_sheet.cal <- data_from_sheet.cal[, -c(6,11)]
data_from_sheet.val <- data_from_sheet.val[, -c(6,11)]

# Data.frame for Histograms
df.flow.OBS.canas <- rbind(data_from_sheet.cal, data_from_sheet.val)
df.flow.OBS.canas <- df.flow.OBS.canas[, -c(1,3,4)]
df.flow.OBS.canas <- melt(df.flow.OBS.canas)

# Calculate exceedance probabilities
data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$observed_Q, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$observed_Q_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_gr2m, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q_gr2m_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

#data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_abcd, decreasing = TRUE), ]
#data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
#data_from_sheet.cal$Q_abcd_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q.budyko, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q.budyko_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_awbm, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q_awbm_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_ihacres, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q_ihacres_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q.simhyd, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q.simhyd_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

#data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$mean, decreasing = TRUE), ]
#data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
#data_from_sheet.cal$mean_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

# Calculate exceedance probabilities
data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$observed_Q, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$observed_Q_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_gr2m, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q_gr2m_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

#data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_abcd, decreasing = TRUE), ]
#data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
#data_from_sheet.val$Q_abcd_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q.budyko, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q.budyko_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_awbm, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q_awbm_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_ihacres, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q_ihacres_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q.simhyd, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q.simhyd_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

#data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$mean, decreasing = TRUE), ]
#data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
#data_from_sheet.val$mean_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

# blended.OBS data.frames
df.blended.OBS <- rbind(data_from_sheet.cal, data_from_sheet.val)
df.blended.OBS <- df.blended.OBS[, -c(1,3,4)]

df.blended.OBS.f1  <- df.blended.OBS[, c(7,8,  1,10)]
df.blended.OBS.f1$model <- c("obs")
names(df.blended.OBS.f1) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f2  <- df.blended.OBS[, c(7,8,  2,11)]
df.blended.OBS.f2$model <- c("GR2M")
names(df.blended.OBS.f2) <- c("catchment", "stage", "Q", "Prob", "model")

#df.blended.OBS.f3  <- df.blended.OBS[, c(9,10,  3,14)]
#df.blended.OBS.f3$model <- c("abcd")
#names(df.blended.OBS.f3) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f4  <- df.blended.OBS[, c(7,8,  3,12)]
df.blended.OBS.f4$model <- c("BUDYKO")
names(df.blended.OBS.f4) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f5  <- df.blended.OBS[, c(7,8,  4,13)]
df.blended.OBS.f5$model <- c("AWBM")
names(df.blended.OBS.f5) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f6  <- df.blended.OBS[, c(7,8,  5,14)]
df.blended.OBS.f6$model <- c("ihaces")
names(df.blended.OBS.f6) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f7  <- df.blended.OBS[, c(7,8,  6,15)]
df.blended.OBS.f7$model <- c("SIMHYD")
names(df.blended.OBS.f7) <- c("catchment", "stage", "Q", "Prob", "model")

#df.blended.OBS.f8  <- df.blended.OBS[, c(9,10,  8,19)]
#df.blended.OBS.f8$model <- c("mean")
#names(df.blended.OBS.f8) <- c("catchment", "stage", "Q", "Prob", "model")

df.melt.abagares <- rbind(df.blended.OBS.f1,
                          df.blended.OBS.f2,
                          #df.blended.OBS.f3,
                          df.blended.OBS.f4,
                          df.blended.OBS.f5,
                          df.blended.OBS.f6,
                          df.blended.OBS.f7)
                          #df.blended.OBS.f8)

#scale02 <- c("#0571b0", "black", "#a6761d","#ca0020","#4dac26", "#2166ac", "blue")
scale02 <- c ("#ffd92f", "#e41a1c","#4daf4a","blue", "black","#ff7f00","#dd1c77", "#ffd92f")

# FDC Flow Duration Curve for canas -------------------------
ggplot() +
  geom_line(aes(x = Prob,y = Q,colour = model,linetype = stage),data=df.melt.abagares,size = 0.95) +
  facet_wrap(facets = ~stage) +
  geom_vline(data=df.melt.abagares,xintercept = 10.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 25.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 50.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 75.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 95.0,colour = '#969696',size = 0.55,linetype = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("canas") +
  xlab("Prob") +
  ylab("Q (mm)") +
  #scale_y_continuous(trans = "log10") +
  #scale_x_continuous(trans = "log10") +
  scale_color_manual(values = scale02) +
  theme_bw()

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Hydrographs
#--------------------------------------------------------------------------------------------------------------

# A quasi-date ID is added
df.subset.length <- subset(df.flow.OBS.canas, variable == "observed_Q")
df.flow.OBS.canas$Serial_Number <- rep((1:nrow(df.subset.length)),6)
df.flow.OBS.canas$test <- rep(((df.subset.length$value)),6)


# Hydrograph test canas -----
ggplot() +
  geom_point(aes(x = Serial_Number,y = value,shape = stage,colour = variable),
             data=df.flow.OBS.canas,size = 0.99) +
  geom_line(aes(x = Serial_Number,y = value,colour = variable,linetype = stage),
            data=df.flow.OBS.canas,size = 0.45) +
  scale_color_manual(values = scale02) +
  ggtitle("canas") +
  xlab("Prob") +
  ylab("Q (mm)") +
  theme_bw()

# Hydrograph test canas -----
ggplot() +
  geom_point(aes(x = Serial_Number,y = value,shape = stage,colour = variable),data=df.flow.OBS.canas,size = 0.99) +
  geom_line(aes(x = Serial_Number,y = test),data=df.flow.OBS.canas) +
  geom_line(aes(x = Serial_Number,y = value,colour = variable,linetype = stage),data=df.flow.OBS.canas,size = 0.75) +
  scale_color_manual(values = scale02) +
  ggtitle("canas") +
  xlab("Prob") +
  ylab("Q (mm)") +
  theme_bw() +
  facet_wrap(facets = ~variable)

scale03 <- c ("#ffd92f", "#e41a1c","#4daf4a","blue", "#dd1c77","#ff7f00","#dd1c77", "#ffd92f")

# Historical DATE vectors are created according to calibration period
dates_obs <- seq(as.Date("1966-01-01"), as.Date("1989-01-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs <- as.POSIXlt(dates_obs, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
df.flow.OBS.canas$DatesR <- dates_obs


# Convert the 'DatesR' column to the Date class
df.flow.OBS.canas$DatesR <- as.Date(df.flow.OBS.canas$DatesR, format = "%Y-%m-%d")

# data.frame is requested
View(df.flow.OBS.canas)

# Variables rename
df.flow.OBS.canas$variable <- ifelse(df.flow.OBS.canas$variable == "Q_gr2m", "GR2M",
                                     ifelse(df.flow.OBS.canas$variable == "Q.budyko", "BUDYKO",
                                            ifelse(df.flow.OBS.canas$variable == "Q_awbm", "AWMB",
                                                   ifelse(df.flow.OBS.canas$variable == "Q.simhyd", "SIMHYD",
                                                          ifelse(df.flow.OBS.canas$variable == "Q_ihacres", "IHACRES",
                                                                 ifelse(df.flow.OBS.canas$variable == "observed_Q", "OBS",1))))))

# Hydrograph test canas -----
ggplot() +
  geom_point(aes(x = DatesR,
                 y = value,
                 shape = stage,
                 colour = stage),
             data=df.flow.OBS.canas,size = 0.99) +
  geom_line(aes(x = DatesR,
                y = value,
                colour = stage,
                linetype = stage),
            data=df.flow.OBS.canas,size = 0.75) +
  geom_line(aes(x = DatesR,
                y = test,
                colour = "orange"),
            data=df.flow.OBS.canas,size = 0.55) +
  # Use scale_x_date() to handle the Date class data
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") + 
  scale_color_manual(values = scale03) +
  scale_colour_manual(values = c("VAL" = "#ca0020", "CAL" = "#2166ac")) +
  scale_fill_manual(values = c("VAL" = "#ca0020", "CAL" = "#2166ac")) +
  theme_bw() +
  ggtitle("Cannas") +
  xlab("Date (Years)") + # Changed label to reflect the X-axis content
  ylab("Streamflow (mm)") +
  facet_grid(rows = variable ~ .)

# An export/import data.frame is created to free RA memory
save(df.flow.OBS.canas, file = "windows.RData")

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: barrancas
#--------------------------------------------------------------------------------------------------------------

#----------- barrancas ----------
workbook <- loadWorkbook("Multi_Hydro_Ensemble_barrancas_SECOND-ROUND.xlsx")
sheet_name.cal <- "df.ensemble.out.cal" # Replace with the name of the sheet you want to export
data_from_sheet.cal <- read.xlsx(workbook, sheet = sheet_name.cal)
sheet_name.val <- "df.ensemble.out.val" # Replace with the name of the sheet you want to export
data_from_sheet.val <- read.xlsx(workbook, sheet = sheet_name.val)
data_from_sheet.cal$catchment <- c("barrancas")
data_from_sheet.val$catchment <- c("barrancas")
data_from_sheet.cal$stage <- c("CAL")
data_from_sheet.val$stage <- c("VAL")
data_from_sheet.val[c(1:47), 13] <- c("CAL")

# abcd and mean models are removed
data_from_sheet.cal <- data_from_sheet.cal[, -c(6,11)]
data_from_sheet.val <- data_from_sheet.val[, -c(6,11)]

# Data.frame for Histograms
df.flow.OBS.barrancas <- rbind(data_from_sheet.cal, data_from_sheet.val)
df.flow.OBS.barrancas <- df.flow.OBS.barrancas[, -c(1,3,4)]
df.flow.OBS.barrancas <- melt(df.flow.OBS.barrancas)

# Calculate exceedance probabilities
data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$observed_Q, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$observed_Q_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_gr2m, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q_gr2m_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

#data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_abcd, decreasing = TRUE), ]
#data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
#data_from_sheet.cal$Q_abcd_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q.budyko, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q.budyko_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_awbm, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q_awbm_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q_ihacres, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q_ihacres_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$Q.simhyd, decreasing = TRUE), ]
data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
data_from_sheet.cal$Q.simhyd_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

#data_from_sheet.cal <- data_from_sheet.cal[order(data_from_sheet.cal$mean, decreasing = TRUE), ]
#data_from_sheet.cal$Rank <- 1:nrow(data_from_sheet.cal)
#data_from_sheet.cal$mean_FDC <- (data_from_sheet.cal$Rank / (nrow(data_from_sheet.cal) + 1)) * 100

# Calculate exceedance probabilities
data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$observed_Q, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$observed_Q_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_gr2m, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q_gr2m_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

#data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_abcd, decreasing = TRUE), ]
#data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
#data_from_sheet.val$Q_abcd_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q.budyko, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q.budyko_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_awbm, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q_awbm_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q_ihacres, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q_ihacres_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$Q.simhyd, decreasing = TRUE), ]
data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
data_from_sheet.val$Q.simhyd_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

#data_from_sheet.val <- data_from_sheet.val[order(data_from_sheet.val$mean, decreasing = TRUE), ]
#data_from_sheet.val$Rank <- 1:nrow(data_from_sheet.val)
#data_from_sheet.val$mean_FDC <- (data_from_sheet.val$Rank / (nrow(data_from_sheet.val) + 1)) * 100

# blended.OBS data.frames
df.blended.OBS <- rbind(data_from_sheet.cal, data_from_sheet.val)
df.blended.OBS <- df.blended.OBS[, -c(1,3,4)]

df.blended.OBS.f1  <- df.blended.OBS[, c(7,8,  1,10)]
df.blended.OBS.f1$model <- c("obs")
names(df.blended.OBS.f1) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f2  <- df.blended.OBS[, c(7,8,  2,11)]
df.blended.OBS.f2$model <- c("GR2M")
names(df.blended.OBS.f2) <- c("catchment", "stage", "Q", "Prob", "model")

#df.blended.OBS.f3  <- df.blended.OBS[, c(9,10,  3,14)]
#df.blended.OBS.f3$model <- c("abcd")
#names(df.blended.OBS.f3) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f4  <- df.blended.OBS[, c(7,8,  3,12)]
df.blended.OBS.f4$model <- c("BUDYKO")
names(df.blended.OBS.f4) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f5  <- df.blended.OBS[, c(7,8,  4,13)]
df.blended.OBS.f5$model <- c("AWBM")
names(df.blended.OBS.f5) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f6  <- df.blended.OBS[, c(7,8,  5,14)]
df.blended.OBS.f6$model <- c("ihaces")
names(df.blended.OBS.f6) <- c("catchment", "stage", "Q", "Prob", "model")

df.blended.OBS.f7  <- df.blended.OBS[, c(7,8,  6,15)]
df.blended.OBS.f7$model <- c("SIMHYD")
names(df.blended.OBS.f7) <- c("catchment", "stage", "Q", "Prob", "model")

#df.blended.OBS.f8  <- df.blended.OBS[, c(9,10,  8,19)]
#df.blended.OBS.f8$model <- c("mean")
#names(df.blended.OBS.f8) <- c("catchment", "stage", "Q", "Prob", "model")

df.melt.abagares <- rbind(df.blended.OBS.f1,
                          df.blended.OBS.f2,
                          #df.blended.OBS.f3,
                          df.blended.OBS.f4,
                          df.blended.OBS.f5,
                          df.blended.OBS.f6,
                          df.blended.OBS.f7)
#df.blended.OBS.f8)

#scale02 <- c("#0571b0", "black", "#a6761d","#ca0020","#4dac26", "#2166ac", "blue")
scale02 <- c ("#ffd92f", "#e41a1c","#4daf4a","blue", "black","#ff7f00","#dd1c77", "#ffd92f")

# FDC Flow Duration Curve for barrancas -------------------------
ggplot() +
  geom_line(aes(x = Prob,y = Q,colour = model,linetype = stage),data=df.melt.abagares,size = 0.95) +
  facet_wrap(facets = ~stage) +
  geom_vline(data=df.melt.abagares,xintercept = 10.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 25.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 50.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 75.0,colour = '#969696',size = 0.55,linetype = 2) +
  geom_vline(data=df.melt.abagares,xintercept = 95.0,colour = '#969696',size = 0.55,linetype = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("barrancas") +
  xlab("Prob") +
  ylab("Q (mm)") +
  #scale_y_continuous(trans = "log10") +
  #scale_x_continuous(trans = "log10") +
  scale_color_manual(values = scale02) +
  theme_bw()

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Hydrographs
#--------------------------------------------------------------------------------------------------------------

# A quasi-date ID is added
df.subset.length <- subset(df.flow.OBS.barrancas, variable == "observed_Q")
df.flow.OBS.barrancas$Serial_Number <- rep((1:nrow(df.subset.length)),6)
df.flow.OBS.barrancas$test <- rep(((df.subset.length$value)),6)


# Hydrograph test barrancas -----
ggplot() +
  geom_point(aes(x = Serial_Number,y = value,shape = stage,colour = variable),
             data=df.flow.OBS.barrancas,size = 0.99) +
  geom_line(aes(x = Serial_Number,y = value,colour = variable,linetype = stage),
            data=df.flow.OBS.barrancas,size = 0.45) +
  scale_color_manual(values = scale02) +
  ggtitle("barrancas") +
  xlab("Prob") +
  ylab("Q (mm)") +
  theme_bw()

# Hydrograph test barrancas -----
ggplot() +
  geom_point(aes(x = Serial_Number,y = value,shape = stage,colour = variable),data=df.flow.OBS.barrancas,size = 0.99) +
  geom_line(aes(x = Serial_Number,y = test),data=df.flow.OBS.barrancas) +
  geom_line(aes(x = Serial_Number,y = value,colour = variable,linetype = stage),data=df.flow.OBS.barrancas,size = 0.75) +
  scale_color_manual(values = scale02) +
  ggtitle("barrancas") +
  xlab("Prob") +
  ylab("Q (mm)") +
  theme_bw() +
  facet_wrap(facets = ~variable)

scale03 <- c ("#ffd92f", "#e41a1c","#4daf4a","blue", "#dd1c77","#ff7f00","#dd1c77", "#ffd92f")


# Historical DATE vectors are created according to calibration period
dates_obs <- seq(as.Date("1977-04-01"), as.Date("1990-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs <- as.POSIXlt(dates_obs, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
df.flow.OBS.barrancas$DatesR <- dates_obs


# Convert the 'DatesR' column to the Date class
df.flow.OBS.barrancas$DatesR <- as.Date(df.flow.OBS.barrancas$DatesR, format = "%Y-%m-%d")

# data.frame is requested
View(df.flow.OBS.barrancas)

# Variables rename
df.flow.OBS.barrancas$variable <- ifelse(df.flow.OBS.barrancas$variable == "Q_gr2m", "GR2M",
                                         ifelse(df.flow.OBS.barrancas$variable == "Q.budyko", "BUDYKO",
                                                ifelse(df.flow.OBS.barrancas$variable == "Q_awbm", "AWMB",
                                                       ifelse(df.flow.OBS.barrancas$variable == "Q.simhyd", "SIMHYD",
                                                              ifelse(df.flow.OBS.barrancas$variable == "Q_ihacres", "IHACRES",
                                                                     ifelse(df.flow.OBS.barrancas$variable == "observed_Q", "OBS",1))))))


# Hydrograph test barrancas -----
ggplot() +
  geom_point(aes(x = DatesR,
                 y = value,
                 shape = stage,
                 colour = stage),
             data=df.flow.OBS.barrancas,size = 0.99) +
  geom_line(aes(x = DatesR,
                y = value,
                colour = stage,
                linetype = stage),
            data=df.flow.OBS.barrancas,size = 0.75) +
  geom_line(aes(x = DatesR,
                y = test,
                colour = "orange"),
            data=df.flow.OBS.barrancas,size = 0.55) +
  # Use scale_x_date() to handle the Date class data
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") + 
  scale_color_manual(values = scale03) +
  scale_colour_manual(values = c("VAL" = "#ca0020", "CAL" = "#2166ac")) +
  scale_fill_manual(values = c("VAL" = "#ca0020", "CAL" = "#2166ac")) +
  theme_bw() +
  ggtitle("Barrancas") +
  xlab("Date (Years)") + # Changed label to reflect the X-axis content
  ylab("Streamflow (mm)") +
  facet_grid(rows = variable ~ .)

# An export/import data.frame is created to free RA memory
save(df.flow.OBS.barrancas, file = "windows.RData")

#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------

































# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Estas tablas son exlusivas para ver la respuesta cuenca por cuenca
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

df.ens.metrics.rankings <- df.ens.metrics.normalization[, c(10,12,13,24)]

#----------- abangares----------
df.subset.abangares.cal <- subset(df.ens.metrics.rankings, catchment == "abangares"  & stage == "CAL")
df.subset.abangares.cal$ranks <- 6-(rank(df.subset.abangares.cal$uni02))
df.subset.abangares.val <- subset(df.ens.metrics.rankings, catchment == "abangares"  & stage == "VAL")
df.subset.abangares.val$ranks <- 6-(rank(df.subset.abangares.val$uni02))
df.abangares.mild <- rbind(df.subset.abangares.cal, df.subset.abangares.val)
#------------------------------

#----------- tenorio----------
df.subset.tenorio.cal <- subset(df.ens.metrics.rankings, catchment == "tenorio"  & stage == "CAL")
df.subset.tenorio.cal$ranks <- 6-(rank(df.subset.tenorio.cal$uni02))
df.subset.tenorio.val <- subset(df.ens.metrics.rankings, catchment == "tenorio"  & stage == "VAL")
df.subset.tenorio.val$ranks <- 6-(rank(df.subset.tenorio.val$uni02))
df.tenorio.mild <- rbind(df.subset.tenorio.cal, df.subset.tenorio.val)
#------------------------------

#----------- tempisque----------
df.subset.tempisque.cal <- subset(df.ens.metrics.rankings, catchment == "tempisque"  & stage == "CAL")
df.subset.tempisque.cal$ranks <- 6-(rank(df.subset.tempisque.cal$uni02))
df.subset.tempisque.val <- subset(df.ens.metrics.rankings, catchment == "tempisque"  & stage == "VAL")
df.subset.tempisque.val$ranks <- 6-(rank(df.subset.tempisque.val$uni02))
df.tempisque.mild <- rbind(df.subset.tempisque.cal, df.subset.tempisque.val)
#------------------------------

#----------- santarosa----------
df.subset.santarosa.cal <- subset(df.ens.metrics.rankings, catchment == "santarosa"  & stage == "CAL")
df.subset.santarosa.cal$ranks <- 6-(rank(df.subset.santarosa.cal$uni02))
df.subset.santarosa.val <- subset(df.ens.metrics.rankings, catchment == "santarosa"  & stage == "VAL")
df.subset.santarosa.val$ranks <- 6-(rank(df.subset.santarosa.val$uni02))
df.santarosa.mild <- rbind(df.subset.santarosa.cal, df.subset.santarosa.val)
#------------------------------

#----------- salitral----------
df.subset.salitral.cal <- subset(df.ens.metrics.rankings, catchment == "salitral"  & stage == "CAL")
df.subset.salitral.cal$ranks <- 6-(rank(df.subset.salitral.cal$uni02))
df.subset.salitral.val <- subset(df.ens.metrics.rankings, catchment == "salitral"  & stage == "VAL")
df.subset.salitral.val$ranks <- 6-(rank(df.subset.salitral.val$uni02))
df.salitral.mild <- rbind(df.subset.salitral.cal, df.subset.salitral.val)
#------------------------------

#----------- morote----------
df.subset.morote.cal <- subset(df.ens.metrics.rankings, catchment == "morote"  & stage == "CAL")
df.subset.morote.cal$ranks <- 6-(rank(df.subset.morote.cal$uni02))
df.subset.morote.val <- subset(df.ens.metrics.rankings, catchment == "morote"  & stage == "VAL")
df.subset.morote.val$ranks <- 6-(rank(df.subset.morote.val$uni02))
df.morote.mild <- rbind(df.subset.morote.cal, df.subset.morote.val)
#------------------------------

#----------- lagarto----------
df.subset.lagarto.cal <- subset(df.ens.metrics.rankings, catchment == "lagarto"  & stage == "CAL")
df.subset.lagarto.cal$ranks <- 6-(rank(df.subset.lagarto.cal$uni02))
df.subset.lagarto.val <- subset(df.ens.metrics.rankings, catchment == "lagarto"  & stage == "VAL")
df.subset.lagarto.val$ranks <- 6-(rank(df.subset.lagarto.val$uni02))
df.lagarto.mild <- rbind(df.subset.lagarto.cal, df.subset.lagarto.val)
#------------------------------

#----------- elsalto----------
df.subset.elsalto.cal <- subset(df.ens.metrics.rankings, catchment == "elsalto"  & stage == "CAL")
df.subset.elsalto.cal$ranks <- 6-(rank(df.subset.elsalto.cal$uni02))
df.subset.elsalto.val <- subset(df.ens.metrics.rankings, catchment == "elsalto"  & stage == "VAL")
df.subset.elsalto.val$ranks <- 6-(rank(df.subset.elsalto.val$uni02))
df.elsalto.mild <- rbind(df.subset.elsalto.cal, df.subset.elsalto.val)
#------------------------------

#----------- corobici----------
df.subset.corobici.cal <- subset(df.ens.metrics.rankings, catchment == "corobici"  & stage == "CAL")
df.subset.corobici.cal$ranks <- 6-(rank(df.subset.corobici.cal$uni02))
df.subset.corobici.val <- subset(df.ens.metrics.rankings, catchment == "corobici"  & stage == "VAL")
df.subset.corobici.val$ranks <- 6-(rank(df.subset.corobici.val$uni02))
df.corobici.mild <- rbind(df.subset.corobici.cal, df.subset.corobici.val)
#------------------------------

#----------- canas ----------
df.subset.canas.cal <- subset(df.ens.metrics.rankings, catchment == "canas"  & stage == "CAL")
df.subset.canas.cal$ranks <- 6-(rank(df.subset.canas.cal$uni02))
df.subset.canas.val <- subset(df.ens.metrics.rankings, catchment == "canas"  & stage == "VAL")
df.subset.canas.val$ranks <- 6-(rank(df.subset.canas.val$uni02))
df.canas.mild <- rbind(df.subset.canas.cal, df.subset.canas.val)
#------------------------------

#----------- blanco----------
df.subset.blanco.cal <- subset(df.ens.metrics.rankings, catchment == "blanco"  & stage == "CAL")
df.subset.blanco.cal$ranks <- 6-(rank(df.subset.blanco.cal$uni02))
df.subset.blanco.val <- subset(df.ens.metrics.rankings, catchment == "blanco"  & stage == "VAL")
df.subset.blanco.val$ranks <- 6-(rank(df.subset.blanco.val$uni02))
df.blanco.mild <- rbind(df.subset.blanco.cal, df.subset.blanco.val)
#------------------------------

#----------- barrancas----------
df.subset.barrancas.cal <- subset(df.ens.metrics.rankings, catchment == "barrancas"  & stage == "CAL")
df.subset.barrancas.cal$ranks <- 6-(rank(df.subset.barrancas.cal$uni02))
df.subset.barrancas.val <- subset(df.ens.metrics.rankings, catchment == "barrancas"  & stage == "VAL")
df.subset.barrancas.val$ranks <- 6-(rank(df.subset.barrancas.val$uni02))
df.barrancas.mild <- rbind(df.subset.barrancas.cal, df.subset.barrancas.val)
#------------------------------

df.calibration.mild <- rbind(df.subset.abangares.cal,
                             df.subset.tenorio.cal,
                             df.subset.tempisque.cal,
                             df.subset.santarosa.cal,
                             df.subset.salitral.cal,
                             df.subset.morote.cal,
                             df.subset.lagarto.cal,
                             df.subset.elsalto.cal,
                             df.subset.corobici.cal,
                             df.subset.canas.cal,
                             df.subset.blanco.cal,
                             df.subset.barrancas.cal)

df.calibration.mild.N1 <- subset(df.calibration.mild, ranks == 1)


df.validation.mild <- rbind(df.subset.abangares.val,
                            df.subset.tenorio.val,
                            df.subset.tempisque.val,
                            df.subset.santarosa.val,
                            df.subset.salitral.val,
                            df.subset.morote.val,
                            df.subset.lagarto.val,
                            df.subset.elsalto.val,
                            df.subset.corobici.val,
                            df.subset.canas.val,
                            df.subset.blanco.val,
                            df.subset.barrancas.val)

df.validation.mild.N1 <- subset(df.validation.mild, ranks == 1)



# An export/import data.frame is created to free RA memory
save(df.ens.metrics , file = "windows.RData")






#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------




