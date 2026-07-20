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

#-------------------------------------------------------------------------------------------------------------------

# Workspace is cleared
gc(); rm(list = ls())

# Scientific notation is disabled
options(scipen=999)

# dredge {MuMIn} configuration is set
options(na.action = "na.fail")

# Start time is recorded
start.time <- Sys.time()

# Working directory is defined
setwd("~/Dropbox/Academics/IDF_CC_tool_CANADA/R_scripts/GR2M_REGIO")

set.seed(1230) # 1237 / 1230

# CRAN libraries are loaded
require(airGR)
require(airGRteaching)
require(DescTools)
require(dplyr)
require(ggplot2)
require(ggpubr)
require(lubridate)
require(matrixStats)
require(pastecs)
require(plyr)
require(reshape)
require(reshape2)
require(tidyr)
require(viridis)
require(weathermetrics)

# Working directory is defined
setwd("~/Dropbox/Academics/IDF_CC_tool_CANADA/R_scripts/GR2M_REGIO")

# Scientific notation is disabled
options(scipen=999)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Watershed corobici
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: MuMin Validation Data
#--------------------------------------------------------------------------------------------------------------

# Regression calibration X1/X2 parameters are loaded
df.mumin <- read.table("mumin_validation_ml_spa_psa.txt",header=T,sep="\t",quote="")

# A catchment subset is executed
df.mumin <- subset(df.mumin, catchment == "corobici")

# Original Michel model are copied to data.frame
df.mumin [12,] <- df.mumin [11,]
df.mumin[12,2] <- df.mumin[12,3]
df.mumin[12,4] <- df.mumin[12,5]
df.mumin[12,6] <- c("original")

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Observed Calibration and Validation
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs <- read.table("airGR_corobici_2025.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs$E <- BasinObs$E*30
BasinObs$P <- BasinObs$P*1.0
BasinObs$Qmm <- BasinObs$Qmm*1.0

# Historical DATE vectors are created according to calibration period
dates_obs <- seq(as.Date("1961-01-01"), as.Date("1990-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs <- as.POSIXlt(dates_obs, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs$DatesR <- dates_obs

# data.frame is requested
View(BasinObs)

# Flowpeaks are compensated due to overestimation
factorC <- 1.00 # 0.0 to 1.00 proportional factor is created

# A temporal vector is created
basin_temp_01 <- BasinObs$Qmm*factorC

# Maximum flowpeak threshold is defined (mm/month)
basin_temp_02 <-ifelse(BasinObs$Qmm >= 999, basin_temp_01, BasinObs$Qmm)

# Compensated flows are replaced
BasinObs$Qmm <- basin_temp_02

# GR2M InputsModel object is created
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                 DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P,
                                 PotEvap = BasinObs$E)

# List containers are erased
list.export <- NULL
list.cont_nse.cal <- NULL
list.cont_kge.cal <- NULL
list.cont_nse.val <- NULL
list.cont_kge.val <- NULL

# List containers are defined
list.cont_nse.cal <- list()
list.cont_kge.cal <- list()
list.cont_nse.val <- list()
list.cont_kge.val <- list()
list.export <- list()

# ==============================
# Outermost loop is initialized
# ==============================

for (cc in 1:12) { # length = 12: months
  
  #--------------------------------------------------------------------------------------------------------------
  # SUBBLOCK: Historical Modelling
  #--------------------------------------------------------------------------------------------------------------
  
  # CALIBRATION ++++++++++
  
  # GR2M run period is selected
  Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m")=="1961-01"),
                 which(format(BasinObs$DatesR, format = "%Y-%m")=="1970-01"))
  
  # GR2M RunOptions object is created
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                 InputsModel = InputsModel,
                                 IndPeriod_Run = Ind_Run)
  
  # Optimum parameters are defined
  Param <- c(X1 =  df.mumin[cc, 2], X2 = df.mumin[cc, 4])
  
  # GR2M OutputsModel object is created
  OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                                RunOptions = RunOptions,
                                Param = Param)
  
  # A GR2M summary plot is requested
  plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run],
       cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
       BasinArea=929.40, which = c("all"))
  
  # Nash-Sutcliffe Efficiency is defined
  InputsCrit  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                  InputsModel = InputsModel,
                                  RunOptions = RunOptions,
                                  Obs = BasinObs$Qmm[Ind_Run])
  
  # Nash-Sutcliffe Efficiency is requested
  OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit,
                               OutputsModel = OutputsModel)
  
  
  # KGE Efficiency is defined
  InputsCrit02  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE,
                                    InputsModel = InputsModel,
                                    RunOptions = RunOptions,
                                    Obs = BasinObs$Qmm[Ind_Run])
  
  # KGE Efficiency is requested
  OutputsCrit02 <- ErrorCrit_KGE(InputsCrit = InputsCrit02,
                                 OutputsModel = OutputsModel)
  
  # VALIDATION ++++++++++
  
  # GR2M run period is selected
  Ind_Run02 <- seq(which(format(BasinObs$DatesR, format = "%Y-%m")=="1970-01"),
                   which(format(BasinObs$DatesR, format = "%Y-%m")=="1979-11"))
  
  # GR2M RunOptions02 object is created
  RunOptions02 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                   InputsModel = InputsModel,
                                   IndPeriod_Run = Ind_Run02)
  
  # GR2M OutputsModel02 object is created
  OutputsModel02 <- RunModel_GR2M(InputsModel = InputsModel,
                                  RunOptions = RunOptions02,
                                  Param = Param)
  
  # A GR2M summary plot is requested
  plot(OutputsModel02, Qobs = BasinObs$Qmm[Ind_Run02],
       cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
       BasinArea=929.40, which = c("all"))
  
  # Nash-Sutcliffe Efficiency is defined
  InputsCritVal  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                     InputsModel = InputsModel,
                                     RunOptions = RunOptions02,
                                     Obs = BasinObs$Qmm[Ind_Run02])
  
  # Nash-Sutcliffe Efficiency is requested
  OutputsCritVal <- ErrorCrit_NSE(InputsCrit = InputsCritVal,
                                  OutputsModel = OutputsModel02)
  
  
  # KGE Efficiency is defined
  InputsCritVal02  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE,
                                       InputsModel = InputsModel,
                                       RunOptions = RunOptions02,
                                       Obs = BasinObs$Qmm[Ind_Run02])
  
  # KGE Efficiency is requested
  OutputsCritVal02 <- ErrorCrit_KGE(InputsCrit = InputsCritVal02,
                                    OutputsModel = OutputsModel02)
  
  
  # NSE metric is defined
  v.nse.cal <- OutputsCrit$CritValue
  v.nse.val <- OutputsCritVal$CritValue

  # KGE metric is defined
  v.kge.cal <- OutputsCrit02$CritValue
  v.kge.val <- OutputsCritVal02$CritValue

  # MD metric is defined
  v.md.cal <- hydroGOF::md(OutputsModel$Qsim, InputsCrit$Obs)
  v.md.val <- hydroGOF::md(OutputsModel02$Qsim, InputsCritVal02$Obs)
  
  # R metric is defined
  v.r.cal <- hydroGOF::rPearson(OutputsModel$Qsim, InputsCrit$Obs)
  v.r.val <- hydroGOF::rPearson(OutputsModel02$Qsim, InputsCritVal02$Obs)
  
  # PBIAS metric is defined
  v.pbias.cal <- hydroGOF::pbias(OutputsModel$Qsim, InputsCrit$Obs)
  v.pbias.val <- hydroGOF::pbias(OutputsModel02$Qsim, InputsCritVal02$Obs)
  
  # VE metric is defined
  v.ve.cal <- hydroGOF::VE(OutputsModel$Qsim, InputsCrit$Obs)
  v.ve.val <- hydroGOF::VE(OutputsModel02$Qsim, InputsCritVal02$Obs)
  
  # LF metric is defined
  v.lf.cal <- hydroGOF::KGElf(OutputsModel$Qsim, InputsCrit$Obs)
  v.lf.val <- hydroGOF::KGElf(OutputsModel02$Qsim, InputsCritVal02$Obs)

  #--------------------------------------------------------------------------------------------------------------
  # SUBBLOCK: Compiling results
  #--------------------------------------------------------------------------------------------------------------
  
  # A list container is defined
  list.export[[cc]] <- as.data.frame(t(c(v.nse.cal,
                                          v.kge.cal,
                                          v.nse.val,
                                          v.kge.val,
                                          v.md.cal,
                                          v.r.cal,
                                          v.md.val,
                                          v.r.val,
                                          v.pbias.cal,
                                          v.pbias.val,
                                          v.ve.cal,
                                          v.ve.val,
                                          v.lf.cal,
                                          v.lf.val,
                                          df.mumin[cc, 1],
                                          df.mumin[cc, 6])))
  
  # ==============================
} # Outermost loop is closed
  # ==============================

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Metrics Development
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# An external list of "date" data.frames is created
df.export <- plyr::ldply(list.export, data.frame)

# data.frame variables are renamed
names(df.export) <- c("nse_cal",
                       "kge_cal",
                       "nse_val",
                       "kge_val",
                       "md_cal",
                       "r_cal",
                       "md_val",
                       "r_val",
                       "pbias_cal",
                       "pbias_val",
                       "ve_cal",
                       "ve_val",
                       "lf_cal",
                       "lf_val",
                       "catchment",
                       "method")

# data.frame variables are converted to numeric
df.export[1:14] <- lapply(df.export[1:14], as.numeric)

# data.frame variables are rounded
df.export[,-c(15,16)] <-round(df.export[,-c(15,16)], 3) #the "-1" excludes column 1

# A data.frame is requested
View(df.export)

# A dotchart is requested
ggdotchart(df.export, y = "kge_cal", x = "method", 
  group = "method", color = "method", palette = "jco", #group = "method",#
  add = "segment", #position = position_dodge(0.3),
  sorting = "descending", #facet.by = "method",
  rotate = TRUE, legend = "none")

#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------

