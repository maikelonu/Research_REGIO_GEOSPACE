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
require(DWBmodelUN) # Concepts behind DWBmodelUN / https://github.com/dazamora/DWBmodelUN
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

# In case these libraries need to be installed:
# devtools::install_github('tanerumit/hydrosystems')
# devtools::install_github("JosephGuillaume/hydromad")
# devtools::install_github("dazamora/DWBmodelUN")

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Weather Station Metadata + Input Data
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Peakflows are compensated due to overestimation
factorC <- 1.0

# Maximum peakflow threshold is defined (mm/month)
factorThres <- 9999

# Watershed compiled files are loaded
BasinObs <- read.table("airGR_tempisque_2025.txt",header=T,sep="\t",quote="")

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

# A temporal vector is created
basin_temp_01 <- BasinObs$Qmm*factorC

# Maximum flowpeak threshold is defined (mm/month)
basin_temp_02 <-ifelse(BasinObs$Qmm > factorThres, basin_temp_01, BasinObs$Qmm)

# Compensated flows are replaced
BasinObs$Qmm <- basin_temp_02

# Calibration run period is selected
Ind.Run.Cal <- seq(which(format(BasinObs$DatesR, format = "%Y-%m")=="1979-10"),
                   which(format(BasinObs$DatesR, format = "%Y-%m")=="1985-03"))

# Validation run period is selected
Ind.Run.Val <- seq(which(format(BasinObs$DatesR, format = "%Y-%m")=="1985-03"),
                   which(format(BasinObs$DatesR, format = "%Y-%m")=="1988-04"))

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: GR2M Hydrological model
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Initial fixed parameters are defined
Param <- c(X1 =  1500, X2 = 0.50)

# Parameter ranges are defined
Param.LOW <- c(X1 =  100, X2 = 0.05)
Param.HIG <- c(X1 =  3000, X2 = 2.0)

# Sensitivity parameter ranges are defined
Param.Ranges <- data.frame(min = c(100, 0.05),
                           max = c(3000, 2.0))

# data.frames variales are defined
rownames(Param.Ranges) <- c("X1", "X2")

# Parameter data containers are defined
pars.INI <- NULL
pars.INI <- list()

# ==============================
# Outermost loop is initialized
# ==============================

# Initial random parameters are defined
for (i in 1:100) {
  set.seed(115 + i)
  pars.INI[i] <- list(c(X1 = runif(1, 100, 3000),
                        X2 = runif(1, 0.05, 2.00)))
  # ============================
} # End of outermost-loop
  # ============================

# GR2M InputsModel object is created
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                 DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P,
                                 PotEvap = BasinObs$E)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Calibration
#--------------------------------------------------------------------------------------------------------------

# Objective function selector is defined
of.selector <- 6 

# Available objective functions:
# sswr =  #3 opti
# nse =   #4 opti
# lnse =  #5 opti
# kge =   #6 opti
# kgenp = #7 opti

# Calibration run period is defined
Ind.Run = Ind.Run.Cal

# GR2M RunOptions object is created
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                               InputsModel = InputsModel,
                               IndPeriod_Run = Ind.Run)

#--------------------------------------------------------------------------------------------------------------
# GR2M hydrological model function is defined for optimization
f.GR2M.opti <- function(Param) {
  
  # GR2M OutputsModel object is created
  OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                                RunOptions = RunOptions,
                                Param = Param)
  
  # Length of observations and simulations is checked
  length(BasinObs$Qmm[Ind.Run])
  length(OutputsModel$Qsim)
  
  # Observations and simulations are isolated
  outQOBS <- BasinObs$Qmm[Ind.Run]
  QSIM <- OutputsModel$Qsim
  
  # Sum of observed runoff [mm]
  of.sum.obs <- sum(outQOBS, na.rm = TRUE)
  #print(paste0("OBS.RUN: ", of.sum.obs))
  
  # Sum of simulated runoff [mm]
  of.sum.sim <- sum(QSIM, na.rm = TRUE)
  #print(paste0("OBS.SIM: ", of.sum.sim))
  
  # SSWR (sum of squared weighted residuals) [mm2] {minimize} (-)
  of.sswr <- sum((outQOBS - QSIM)^2, na.rm = TRUE)
  #print(paste0("SSWR: ", of.sswr))
  
  # NSE (Nash-Sutcliffe Efficiency) [-] {maximize} (0-1)
  of.nse <- (1-(hydroGOF::NSE(QSIM, outQOBS)))
  #print(paste0("NSE: ", of.nse))
  
  # LNSE (Log of the Nash-Sutcliffe Efficiency) [-] {maximize} (0-1)
  of.lnse <- (abs(hydroGOF::NSE(log(QSIM), log(outQOBS))))
  #print(paste0("LNSE: ", of.lnse))
  
  # KGE (Kling-Gupta Efficiency) [-] {maximize} (0-1)
  of.kge <- (1-(hydroGOF::KGE(QSIM, outQOBS)))
  #print(paste0("KGE: ", of.kge))
  
  # KGEnp (Non-parametric version of the Kling-Gupta Efficiency) [-] {maximize} (0-1)
  of.kgenp <- (1-(hydroGOF::KGEnp(QSIM, outQOBS)))
  #print(paste0("KGEnp: ", of.kgenp ))
  
  # PBIAS (Percent Bias) {minimize} (-)
  of.pbias <- (hydroGOF::pbias(outQOBS, QSIM, na.rm=TRUE))
  #print(paste0("PBIAS: ", of.pbias))
  
  # IA (Index_Agreement) {maximize} (0-1)
  of.ia <- (hydroGOF::d(outQOBS, QSIM))
  #print(paste0("IA: ", of.ia))
  
  # MIA (Modified_Index_Agreement) {maximize} (0-1)
  of.mia <- (hydroGOF::md(outQOBS, QSIM))
  #print(paste0("MIA: ", of.mia))
  
  # RIA (Refined_Index_Agreement) {maximize} (0-1)
  of.rsr <- (hydroGOF::rsr(outQOBS, QSIM))
  #print(paste0("RIA: ", of.rsr))
  
  # Pearson Correlation Coefficient (RPearson) {maximize} (0-1)
  of.corr <- (hydroGOF::R2(outQOBS, QSIM))
  #print(paste0("RPearson: ", of.corr))
  
  # Root Mean Square Error (RMSE) {minimize} (-)
  of.rmse <- (rmse(outQOBS, QSIM))
  #print(paste0("RMSE: ", of.rmse))
  
  # Multi Objective Functions
  of.test.01 <- 0.05*(of.nse) + 0.95*(of.nse)
  #print(paste0("of.test.01: ",of.test.01))
  
  # Output vector is created
  vector.sum <- c(of.sum.obs, of.sum.sim, of.sswr, of.nse, 
                  of.lnse, of.kge, of.kgenp, of.pbias,
                  of.ia, of.mia, of.rsr, of.corr, of.rmse, of.test.01)
  
  # Output data.frame is created
  df.sum <- data.frame(sum_obs = c(of.sum.obs),      #1
                                  sum_sim = c(of.sum.sim),      #2
                                  sswr = c(of.sswr),            #3 opti
                                  nse = c(of.nse),              #4 opti
                                  lnse = c(of.lnse),            #5 opti
                                  kge = c(of.kge),              #6 opti
                                  kgenp = c(of.kgenp),          #7 opti
                                  pbias = c(of.pbias),          #8
                                  ia = c(of.ia),                #9
                                  mia = c(of.mia),              #10
                                  rsr = c(of.rsr),              #11
                                  corr = c(of.corr),            #12
                                  rmse = c(of.rmse),            #13
                                  of_test_01 = c(of.test.01))   #14
  
  # Function return is defined
  return(as.vector(df.sum[1, of.selector]))
}
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# GR2M hydrological model function is defined for production
f.GR2M.prod <- function(Param) {
  
  # GR2M OutputsModel object is created
  OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                                RunOptions = RunOptions,
                                Param = Param)
  
  # Length of observations and simulations is checked
  length(BasinObs$Qmm[Ind.Run])
  length(OutputsModel$Qsim)
  
  # Observations and simulations are isolated
  outQOBS <- BasinObs$Qmm[Ind.Run]
  QSIM <- OutputsModel$Qsim
  
  # Sum of observed runoff [mm]
  of.sum.obs <- sum(outQOBS, na.rm = TRUE)
  #print(paste0("OBS.RUN: ", of.sum.obs))
  
  # Sum of simulated runoff [mm]
  of.sum.sim <- sum(QSIM, na.rm = TRUE)
  #print(paste0("OBS.SIM: ", of.sum.sim))
  
  # SSWR (sum of squared weighted residuals) [mm2] {minimize} (-)
  of.sswr <- sum((outQOBS - QSIM)^2, na.rm = TRUE)
  #print(paste0("SSWR: ", of.sswr))
  
  # NSE (Nash-Sutcliffe Efficiency) [-] {maximize} (0-1)
  of.nse <- (1-(hydroGOF::NSE(QSIM, outQOBS)))
  #print(paste0("NSE: ", of.nse))
  
  # LNSE (Log of the Nash-Sutcliffe Efficiency) [-] {maximize} (0-1)
  of.lnse <- (abs(hydroGOF::NSE(log(QSIM), log(outQOBS))))
  #print(paste0("LNSE: ", of.lnse))
  
  # KGE (Kling-Gupta Efficiency) [-] {maximize} (0-1)
  of.kge <- (1-(hydroGOF::KGE(QSIM, outQOBS)))
  #print(paste0("KGE: ", of.kge))
  
  # KGEnp (Non-parametric version of the Kling-Gupta Efficiency) [-] {maximize} (0-1)
  of.kgenp <- (1-(hydroGOF::KGEnp(QSIM, outQOBS)))
  #print(paste0("KGEnp: ", of.kgenp ))
  
  # PBIAS (Percent Bias) {minimize} (-)
  of.pbias <- (hydroGOF::pbias(outQOBS, QSIM, na.rm=TRUE))
  #print(paste0("PBIAS: ", of.pbias))
  
  # IA (Index_Agreement) {maximize} (0-1)
  of.ia <- (hydroGOF::d(outQOBS, QSIM))
  #print(paste0("IA: ", of.ia))
  
  # MIA (Modified_Index_Agreement) {maximize} (0-1)
  of.mia <- (hydroGOF::md(outQOBS, QSIM))
  #print(paste0("MIA: ", of.mia))
  
  # RIA (Refined_Index_Agreement) {maximize} (0-1)
  of.rsr <- (hydroGOF::rsr(outQOBS, QSIM))
  #print(paste0("RIA: ", of.rsr))
  
  # Pearson Correlation Coefficient (RPearson) {maximize} (0-1)
  of.corr <- (hydroGOF::R2(outQOBS, QSIM))
  #print(paste0("RPearson: ", of.corr))
  
  # Root Mean Square Error (RMSE) {minimize} (-)
  of.rmse <- (rmse(outQOBS, QSIM))
  #print(paste0("RMSE: ", of.rmse))
  
  # Multi Objective Functions
  of.test.01 <- 0.05*(of.nse) + 0.95*(of.nse)
  #print(paste0("of.test.01: ",of.test.01))
  
  # Output vector is created
  vector.sum <- c(of.sum.obs, of.sum.sim, of.sswr, (1- of.nse), 
                  (1- of.lnse), (1-of.kge), (1- of.kgenp), of.pbias,
                  of.ia, of.mia, of.rsr, of.corr, of.rmse, of.test.01)
  
  # Output data.frame is created
  df.sum <- data.frame(sum_obs = c(of.sum.obs),      #1
                       sum_sim = c(of.sum.sim),      #2
                       sswr = c(of.sswr),            #3 opti
                       nse = c(1 - of.nse),          #4 opti
                       lnse = c(1 - of.lnse),        #5 opti
                       kge = c(1 - of.kge),          #6 opti
                       kgenp = c(1 - of.kgenp),      #7 opti
                       pbias = c(of.pbias),          #8
                       ia = c(of.ia),                #9
                       mia = c(of.mia),              #10
                       rsr = c(of.rsr),              #11
                       corr = c(of.corr),            #12
                       rmse = c(of.rmse),            #13
                       of_test_01 = c(of.test.01))   #14

  # Function return is defined
  return(df.sum)
}
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# GR2M hydrological model function is defined for sensitivity analysis
f.GR2M.sens <- function(Param) {
  
  # GR2M OutputsModel object is created
  OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                                RunOptions = RunOptions,
                                Param = Param)
  
  # Length of observations and simulations is checked
  length(BasinObs$Qmm[Ind.Run])
  length(OutputsModel$Qsim)
  
  # Observations and simulations are isolated
  outQOBS <- BasinObs$Qmm[Ind.Run]
  QSIM <- OutputsModel$Qsim
  
  # Sum of observed runoff [mm]
  of.sum.obs <- sum(outQOBS, na.rm=TRUE)
  #print(paste0("OBS.RUN: ", of.sum.obs))
  
  # Sum of simulated runoff [mm]
  of.sum.sim <- sum(QSIM, na.rm=TRUE)
  #print(paste0("OBS.SIM: ", of.sum.sim))
  
  # SSWR (sum of squared weighted residuals) [mm2] {minimize} (-)
  of.sswr <- sum((outQOBS - QSIM)^2, na.rm = TRUE)
  #print(paste0("SSWR: ", of.sswr))
  
  # NSE (Nash-Sutcliffe Efficiency) [-] {maximize} (0-1)
  of.nse <- (1-(hydroGOF::NSE(QSIM, outQOBS)))
  #print(paste0("NSE: ", of.nse))
  
  # LNSE (Log of the Nash-Sutcliffe Efficiency) [-] {maximize} (0-1)
  of.lnse <- (abs(hydroGOF::NSE(log(QSIM), log(outQOBS))))
  #print(paste0("LNSE: ", of.lnse))
  
  # KGE (Kling-Gupta Efficiency) [-] {maximize} (0-1)
  of.kge <- (1-(hydroGOF::KGE(QSIM, outQOBS)))
  #print(paste0("KGE: ", of.kge))
  
  # KGEnp (Non-parametric version of the Kling-Gupta Efficiency) [-] {maximize} (0-1)
  of.kgenp <- (1-(hydroGOF::KGEnp(QSIM, outQOBS)))
  #print(paste0("KGEnp: ", of.kgenp ))
  
  # PBIAS (Percent Bias) {minimize} (-)
  of.pbias <- (hydroGOF::pbias(outQOBS, QSIM, na.rm=TRUE))
  #print(paste0("PBIAS: ", of.pbias))
  
  # IA (Index_Agreement) {maximize} (0-1)
  of.ia <- (hydroGOF::d(outQOBS, QSIM))
  #print(paste0("IA: ", of.ia))
  
  # MIA (Modified_Index_Agreement) {maximize} (0-1)
  of.mia <- (hydroGOF::md(outQOBS, QSIM))
  #print(paste0("MIA: ", of.mia))
  
  # RIA (Refined_Index_Agreement) {maximize} (0-1)
  of.rsr <- (hydroGOF::rsr(outQOBS, QSIM))
  #print(paste0("RIA: ", of.rsr))
  
  # Pearson Correlation Coefficient (RPearson) {maximize} (0-1)
  of.corr <- (hydroGOF::R2(outQOBS, QSIM))
  #print(paste0("RPearson: ", of.corr))
  
  # Root Mean Square Error (RMSE) {minimize} (-)
  of.rmse <- (rmse(outQOBS, QSIM))
  #print(paste0("RMSE: ", of.rmse))
  
  # Multi Objective Functions
  of.test.01 <- 0.05*(of.nse) + 0.95*(of.nse)
  #print(paste0("of.test.01: ",of.test.01))
  
  # Output vector is created
  vector.sum <- c(of.sum.obs, of.sum.sim, of.sswr, (1- of.nse), 
                  (1- of.lnse), (1-of.kge), (1- of.kgenp), of.pbias,
                  of.ia, of.mia, of.rsr, of.corr, of.rmse, of.test.01)
  
  # Output data.frame is created
  df.sum <- data.frame(sum_obs = c(of.sum.obs),      #1
                       sum_sim = c(of.sum.sim),      #2
                       sswr = c(of.sswr),            #3 opti
                       nse = c(1 - of.nse),          #4 opti
                       lnse = c(1 - of.lnse),        #5 opti
                       kge = c(1 - of.kge),          #6 opti
                       kgenp = c(1 - of.kgenp),      #7 opti
                       pbias = c(of.pbias),          #8
                       ia = c(of.ia),                #9
                       mia = c(of.mia),              #10
                       rsr = c(of.rsr),              #11
                       corr = c(of.corr),            #12
                       rmse = c(of.rmse),            #13
                       of_test_01 = c(of.test.01))   #14

  # Function return is defined
  return(df.sum)
}

#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# Optimization via {optim} quasi-Newton method variable metric algorithm ()
#--------------------------------------------------------------------------------------------------------------

# Sequence of initial chain of parameter set is selected
loop.i = 25

# Optimization object is created
CNT <- optim(par = unlist(pars.INI[loop.i]),
             fn = f.GR2M.opti, 
             method ="L-BFGS-B",
             lower = Param.LOW,
             upper = Param.HIG,
             control = list(maxit = 500))

# Optimization results are requested
CNT$par

# GR2M hydrological model is run with optimum parameters
df.GR2M.out.optim <- f.GR2M.prod(CNT$par)

# Final OFs are requested
vector.gr2m.optim.kge <- df.GR2M.out.optim$kge
vector.gr2m.optim.nse <- df.GR2M.out.optim$nse
vector.gr2m.optim.pbias <- df.GR2M.out.optim$pbias
vector.gr2m.optim.corr <- df.GR2M.out.optim$corr
vector.gr2m.optim.d <- df.GR2M.out.optim$ia
vector.gr2m.optim.md <- df.GR2M.out.optim$mia

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- df.GR2M.out.optim$sswr
N <- length(Ind.Run.Cal)
K <- 2
vector.gr2m.optim.aic <- AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
vector.gr2m.optim.bic <- BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
vector.gr2m.optim.rsr <- df.GR2M.out.optim$rsr
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.gr2m.optim <- as.data.frame(t(c(vector.gr2m.optim.kge,
                                       vector.gr2m.optim.nse,
                                       vector.gr2m.optim.pbias,
                                       vector.gr2m.optim.corr,
                                       vector.gr2m.optim.d,
                                       vector.gr2m.optim.md,
                                       vector.gr2m.optim.aic,
                                       vector.gr2m.optim.bic,
                                       vector.gr2m.optim.rsr,
                                       c("gr2m"),
                                       c("optim"),
                                       c("cal"))))

# data.frame variables are renamed
names(df.cont.gr2m.optim) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

#--------------------------------------------------------------------------------------------------------------
# Optimization via GenSA (Generalized Simulated Annealing Function)
#--------------------------------------------------------------------------------------------------------------

# Sequence of initial chain of parameter set is selected
loop.i = 25

# Optimization object is created
CNT.gr2m <- GenSA(par = unlist(pars.INI[loop.i]),
                  fn = f.GR2M.opti, 
                  lower = Param.LOW,
                  upper = Param.HIG,
                  control = list(maxit = 1000,
                                 nb.stop.improvement = 500,
                                 max.call = 20000,
                                 verbose = TRUE))

# Optimization results are requested
CNT.gr2m$par
CNT.gr2m$counts

# Final optimum parameters data.frame is created
df.gr2m.final.par <- data.frame(as.list(CNT.gr2m$par), as.list(CNT.gr2m$counts)) 

# GR2M hydrological model is run with optimum parameters
df.GR2M.out.gensa <- f.GR2M.prod(CNT.gr2m$par)

# Optimum parameters are converted to data.frame
df.par.out <- (as.data.frame(CNT.gr2m[c("par")]))

# data.frames are rbinded
df.GR2M.out.gensa$X1 <- df.par.out[1,1]
df.GR2M.out.gensa$X2 <- df.par.out[2,1]

# data.frame variables are renamed
df.GR2M.out.gensa$selector <- as.character(of.selector)

# A final optimization test is performed
Param.test.gensa <- c(X1 =  df.GR2M.out.gensa$X1, X2 = df.GR2M.out.gensa$X2)
output.GR2M.opti <- f.GR2M.opti(Param.test.gensa)

# Final production GR2M simulation is executed using Production-Function
Param.prod.gensa <- c(X1 =  df.GR2M.out.gensa$X1, X2 = df.GR2M.out.gensa$X2)
df.GR2M.out.gensa <- f.GR2M.prod(Param.test.gensa)

# Final OFs are requested
vector.gr2m.gensa.kge <- df.GR2M.out.gensa$kge
vector.gr2m.gensa.nse <- df.GR2M.out.gensa$nse
vector.gr2m.gensa.pbias <- df.GR2M.out.gensa$pbias
vector.gr2m.gensa.corr <- df.GR2M.out.gensa$corr
vector.gr2m.gensa.d <- df.GR2M.out.gensa$ia
vector.gr2m.gensa.md <- df.GR2M.out.gensa$mia

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- df.GR2M.out.gensa$sswr
N <- length(Ind.Run.Cal)
K <- 2
vector.gr2m.gensa.aic <- AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
vector.gr2m.gensa.bic <- BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
vector.gr2m.gensa.rsr <- df.GR2M.out.gensa$rsr
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.gr2m.gensa <- as.data.frame(t(c(vector.gr2m.gensa.kge,
                                        vector.gr2m.gensa.nse,
                                        vector.gr2m.gensa.pbias,
                                        vector.gr2m.gensa.corr,
                                        vector.gr2m.gensa.d,
                                        vector.gr2m.gensa.md,
                                        vector.gr2m.gensa.aic,
                                        vector.gr2m.gensa.bic,
                                        vector.gr2m.gensa.rsr,
                                        c("gr2m"),
                                        c("gensa"),
                                        c("cal"))))

# data.frame variables are renamed
names(df.cont.gr2m.gensa) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic","rsr", "model", "method", "stage")

# Optimization results data.frame is created
df.gr2m.final.ofs <- rbind(df.cont.gr2m.optim, df.cont.gr2m.gensa)

# data.frame variables are converted to numeric
df.gr2m.final.ofs[1:9] <- lapply(df.gr2m.final.ofs[1:9], as.numeric)

# Optimization results data.frame is created
t1 <- data.frame(t(c(CNT$par, "gr2m", "optim" )))
names(t1) <- c("x1", "x2", "model", "method")
t2 <- data.frame(t(c(df.par.out[1,1], df.par.out[2,1], "gr2m", "gensa")))
names(t2) <- c("x1", "x2", "model", "method")
df.gr2m.final.par <- rbind(t1, t2)

# data.frame variables are converted to numeric
df.gr2m.final.par[1:2] <- lapply(df.gr2m.final.par[1:2], as.numeric)

# Optimization model calls are added to data.frame
df.gr2m.final.par$call <- as.numeric(CNT.gr2m$counts)

# data.frames are rbinded
X1.final <- df.par.out[1,1]
X2.final <- df.par.out[2,1]

# GR2M hydrological model function is defined for generation
f.GR2M.gen <- function(Param) {
  
  # GR2M OutputsModel object is created
  OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                                RunOptions = RunOptions,
                                Param = Param)
  
  # Length of observations and simulations is checked
  length(BasinObs$Qmm[Ind.Run])
  length(OutputsModel$Qsim)
  
  # Observations and simulations are isolated
  outQOBS <- BasinObs$Qmm[Ind.Run]
  QSIM <- OutputsModel$Qsim
  
  
  # Output data.frame is created
  df.gen <- data.frame(outQOBS, QSIM)   #14
  
  # Function return is defined
  return(df.gen)
}

# GR2M hydrological model function is executed for generation
df.GR2M.gen <- f.GR2M.gen(c(X1.final, X2.final))

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# Calibration run period is defined
Ind.Run = Ind.Run.Val

# GR2M RunOptions object is created
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                               InputsModel = InputsModel,
                               IndPeriod_Run = Ind.Run)

# Final production GR2M simulation is executed using Production-Function
Param.prod.gensa <- c(X1 =  df.gr2m.final.par$x1, X2 = df.gr2m.final.par$x2)
df.GR2M.out.gensa <- f.GR2M.prod(Param.test.gensa)

# Final OFs are requested
vector.gr2m.gensa.kge <- df.GR2M.out.gensa$kge
vector.gr2m.gensa.nse <- df.GR2M.out.gensa$nse
vector.gr2m.gensa.pbias <- df.GR2M.out.gensa$pbias
vector.gr2m.gensa.corr <- df.GR2M.out.gensa$corr
vector.gr2m.gensa.d <- df.GR2M.out.gensa$ia
vector.gr2m.gensa.md <- df.GR2M.out.gensa$mia

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- df.GR2M.out.gensa$sswr
N <- length(Ind.Run.Val)
K <- 2
vector.gr2m.gensa.aic <- AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
vector.gr2m.gensa.bic <- BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
vector.gr2m.gensa.rsr <- df.GR2M.out.gensa$rsr
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.gr2m.gensa.val <- as.data.frame(t(c(vector.gr2m.gensa.kge,
                                            vector.gr2m.gensa.nse,
                                            vector.gr2m.gensa.pbias,
                                            vector.gr2m.gensa.corr,
                                            vector.gr2m.gensa.d,
                                            vector.gr2m.gensa.md,
                                            vector.gr2m.gensa.aic,
                                            vector.gr2m.gensa.bic,
                                            vector.gr2m.gensa.rsr,
                                            c("gr2m"),
                                            c("gensa"),
                                            c("val"))))

# data.frame variables are renamed
names(df.cont.gr2m.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# data.frame variables are converted to numeric
df.cont.gr2m.gensa.val[1:9] <- lapply(df.cont.gr2m.gensa.val[1:9], as.numeric)

# Optimization results data.frame is created
df.gr2m.final.ofs.def <- rbind(df.gr2m.final.ofs, df.cont.gr2m.gensa.val)

# GR2M hydrological model function is executed for generation
df.GR2M.gen.val <- f.GR2M.gen(Param.test.gensa)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: ABCD Hydrological model
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# SOURCES:
# https://abcd.walkerenvres.com/theory.html

# Calibration run period is defined
Ind.Run = Ind.Run.Cal

# Initial fixed parameters are defined
parm.INI <- c(a =  0.950, b = 350.00, c = 0.300,  d = 1.000, S_ini = 1000.00, G_ini = 60.00)

# Parameter ranges are defined
parm.LOW <- c(a =  0.0001, b = 10.000, c = 0.0001,  d = 0.0001, S_ini = 10.00, G_ini = 10.00)
parm.HIG <- c(a =  1.0000, b = 2000.0, c = 1.0000,  d = 1.0000, S_ini = 2000.0, G_ini = 100.00)

#--------------------------------------------------------------------------------------------------------------
# ABCD hydrological model function is defined for production
f.ABCD <- function(parm) {
  
  # ABCD input objects are created
  P.abcd <- BasinObs$P[Ind.Run]
  PE.abcd <- BasinObs$E[Ind.Run]
  Qobs.abcd <- BasinObs$Qmm[Ind.Run]
  
  # abcdQest {hydrosystems} function is requested
  acbd.out <- abcdQest(parm = parm[1:4],
                       P = P.abcd,
                       PE = PE.abcd,
                       S_ini = parm[5],
                       G_ini = parm[6],
                       print.all = FALSE)
  
  # ABCD output objects are created
  acbd.out <- as.vector(acbd.out)
  
  # Metric functions are defined
  of.sse <- sum((Qobs.abcd - acbd.out)^2,     na.rm = TRUE)
  of.nse <- (1-(hydroGOF::NSE(acbd.out, Qobs.abcd)))
  of.kge <- (1-(hydroGOF::KGE(acbd.out, Qobs.abcd)))
  of.kgenp <- (1-(hydroGOF::KGEnp(acbd.out, Qobs.abcd)))
  of.lnse <- (abs(hydroGOF::NSE(log(acbd.out), log(Qobs.abcd))))
  of.test01 <- 0.05*(of.nse) + 0.95*(of.nse)
  
  # sumQ <- sum(acbd.out)
  # df.sum <- data.frame(c(sum(acbd.out)/sum(Qobs.abcd)))
  # names(df.sum) <- c("abcd_out")
  # return(df.sum)
  
  # Function return is defined
  return(of.kge)
}

# An initial simulation is executed
f.ABCD(parm = parm.INI)

# Optimization by GenSA
CNT.abcd <- GenSA(par = parm.INI,
                  fn = f.ABCD, 
                  lower = parm.LOW,
                  upper = parm.HIG,
                  control = list(maxit = 1000,
                                 nb.stop.improvement = 500,
                                 max.call = 20000,
                                 verbose = TRUE))

# Optimum ABCD parameters are requested
CNT.abcd$par
CNT.abcd$counts

# Final optimum parameters data.frame is created
df.abcd.final.par <- data.frame(as.list(CNT.abcd$par), as.list(CNT.abcd$counts)) 

# ABCD model is rerun using optimal parameters
parm.exe.abcd <- c(a = CNT.abcd$par[1],   #1
                   b = CNT.abcd$par[2],   #2
                   c = CNT.abcd$par[3],   #3
                   d = CNT.abcd$par[4],   #4
                   S_ini = CNT.abcd$par[5],
                   G_ini = CNT.abcd$par[6])

# Final KGE is requested
1 - f.ABCD(parm = parm.exe.abcd)

# Optimal parameters are requested
print(parm.exe.abcd)

# A final ABCD run is executed
acbd.out <- abcdQest(parm = parm.exe.abcd[1:4],
                     P = BasinObs$P[Ind.Run.Cal],
                     PE = BasinObs$E[Ind.Run.Cal],
                     S_ini = parm.exe.abcd[5],
                     G_ini = parm.exe.abcd[6],
                     print.all = FALSE)

# Model outputs are isolated
acbd.out <- as.vector(acbd.out)

# Final KGE/NSE are requested
hydroGOF::KGE(acbd.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::NSE(acbd.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::R2(acbd.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::pbias(acbd.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::rPearson(acbd.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::d(acbd.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::md(acbd.out, BasinObs$Qmm[Ind.Run.Cal])

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((acbd.out - BasinObs$Qmm[Ind.Run.Cal])^2, na.rm = TRUE)
N <- length(acbd.out)
K <- 6
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.abcd.gensa.cal <- as.data.frame(t(c(hydroGOF::KGE(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            hydroGOF::NSE(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            hydroGOF::pbias(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            hydroGOF::rPearson(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            hydroGOF::d(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            hydroGOF::md(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            AIC.i,
                                            BIC.i,
                                            hydroGOF::rsr(acbd.out, BasinObs$Qmm[Ind.Run.Cal]),
                                            c("abcd"),
                                            c("gensa"),
                                            c("cal"))))

# data.frame variables are renamed
names(df.cont.abcd.gensa.cal) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# A final ABCD run is executed
acbd.out.val <- abcdQest(parm = parm.exe.abcd[1:4],
                         P = BasinObs$P[Ind.Run.Val],
                         PE = BasinObs$E[Ind.Run.Val],
                         S_ini = parm.exe.abcd[5],
                         G_ini = parm.exe.abcd[6],
                         print.all = FALSE)

# Model outputs are isolated
acbd.out.val <- as.vector(acbd.out.val)

# Final KGE/NSE are requested
hydroGOF::KGE(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::NSE(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::R2(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::pbias(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::rPearson(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::d(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::md(acbd.out.val, BasinObs$Qmm[Ind.Run.Val])

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((acbd.out.val - BasinObs$Qmm[Ind.Run.Val])^2, na.rm = TRUE)
N <- length(acbd.out.val)
K <- 6
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.abcd.gensa.val <- as.data.frame(t(c(hydroGOF::KGE(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            hydroGOF::NSE(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            hydroGOF::pbias(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            hydroGOF::rPearson(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            hydroGOF::d(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            hydroGOF::md(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            AIC.i,
                                            BIC.i,
                                            hydroGOF::rsr(acbd.out.val, BasinObs$Qmm[Ind.Run.Val]),
                                            c("abcd"),
                                            c("gensa"),
                                            c("val"))))

# data.frame variables are renamed
names(df.cont.abcd.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Optimization results data.frame is created
df.abcd.final.ofs.def <- rbind(df.cont.abcd.gensa.cal, df.cont.abcd.gensa.val)

# data.frame variables are converted to numeric
df.abcd.final.ofs.def[1:9] <- lapply(df.abcd.final.ofs.def[1:9], as.numeric)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Budyko Framework Hydrological Model (budyko)
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Calibration run period is defined
Ind.Run = Ind.Run.Cal

# Budyko local function is defined
f.budyko.local <- function (p_v, pet_v, g_v, s_v, alpha1_v, alpha2_v, smax_v, d_v) {
  ncells <- NROW(p_v)
  nmonths <- NCOL(p_v)
  M.result <- matrix(nrow = ncells, ncol = nmonths)
  xo <- M.result
  x <- M.result
  qd <- M.result
  w <- M.result
  yo <- M.result
  y <- M.result
  r <- M.result
  aet <- M.result
  s <- M.result
  qb <- M.result
  g <- M.result
  q_total <- M.result
  dummy <- smax_v - s_v
  xo[, 1] <- dummy + pet_v[, 1]
  x[, 1] <- p_v[, 1] * funFU(PET = xo[, 1], P = p_v[, 1], 
                             alpha = alpha1_v)
  qd[, 1] <- p_v[, 1] - x[, 1]
  w[, 1] <- x[, 1] + s_v
  yo[, 1] <- pet_v[, 1] + smax_v
  y[, 1] <- w[, 1] * funFU(PET = yo[, 1], P = w[, 1], alpha = alpha2_v)
  r[, 1] <- w[, 1] - y[, 1]
  aet[, 1] <- w[, 1] * funFU(PET = pet_v[, 1], P = w[, 
                                                     1], alpha = alpha2_v)
  s[, 1] <- y[, 1] - aet[, 1]
  qb[, 1] <- d_v * g_v
  g[, 1] <- (1 - d_v) * g_v + r[, 1]
  q_total[, 1] <- qb[, 1] + qd[, 1]
  #pb <- txtProgressBar(min = 0, max = nmonths, style = 3)
  for (i in 2:nmonths) {
    dummy <- smax_v - s[, (i - 1)]
    xo[, i] <- dummy + pet_v[, i]
    x[, i] <- p_v[, i] * funFU(PET = xo[, i], P = p_v[, 
                                                      i], alpha = alpha1_v)
    qd[, i] <- p_v[, i] - x[, i]
    w[, i] <- x[, i] + s[, (i - 1)]
    yo[, i] <- pet_v[, i] + smax_v
    y[, i] <- w[, i] * funFU(PET = yo[, i], P = w[, i], 
                             alpha = alpha2_v)
    r[, i] <- w[, i] - y[, i]
    aet[, i] <- w[, i] * funFU(PET = pet_v[, i], P = w[, 
                                                       i], alpha = alpha2_v)
    s[, i] <- y[, i] - aet[, i]
    qb[, i] <- d_v * g[, (i - 1)]
    g[, i] <- (1 - d_v) * g_v + r[, i]
    q_total[, i] <- qb[, i] + qd[, i]
    #setTxtProgressBar(pb, i)
  }
  #close(pb)
  dwb_aux <- list(q_total = q_total, aet = aet, r = r, qd = qd, 
                  qb = qb, s = s, g = g)
  return(dwb_aux)
}

# Initial fixed parameters are defined
parm.bud.INI <- c(g_v =  60, s_v  = 60, alpha1_v  = 0.25,  alpha2_v  = 0.25, smax_v  = 500, d_v = 0.5)

# Parameter ranges are defined
parm.bud.LOW <- c(g_v =  0.00010, s_v  = 0.00010, alpha1_v  = 0.00010,  alpha2_v  = 0.00010, smax_v  = 0.00010, d_v = 0.00001)
parm.bud.HIG <- c(g_v =  100, s_v  = 100, alpha1_v  = 10,  alpha2_v  = 10, smax_v  = 2000.0, d_v = 1)

#--------------------------------------------------------------------------------------------------------------
# budyko hydrological model function is defined for production
f.budyko <- function(parm) {

  # Budyko input objects are created
  P.budyko <- c(BasinObs$P[Ind.Run])
  PE.budyko <- c(BasinObs$E[Ind.Run])
  prueba.P <- data.frame(as.list(P.budyko))
  prueba.PET <- data.frame(as.list(PE.budyko))
  prueba.P <- unname(prueba.P)
  prueba.PET <- unname(prueba.PET)

  # f.budyko.local {DWBmodelUN} function is requested
  budyko.out <-  f.budyko.local(prueba.P, 
                                prueba.PET,
                                g_v = parm[1],
                                s_v = parm[2],
                                alpha1_v = parm[3],
                                alpha2_v = parm[4],
                                smax_v = parm[5],
                                d_v = parm[6])
  
  # budyko output objects are created
  Qo.budyko.vector <- as.vector(budyko.out$q_total)

  # Metric functions are defined
  of.sse <- sum((Qo.budyko.vector - BasinObs$Qmm[Ind.Run])^2,     na.rm = TRUE)
  of.nse <- (1-(hydroGOF::NSE(Qo.budyko.vector, BasinObs$Qmm[Ind.Run])))
  of.kge <- (1-(hydroGOF::KGE(Qo.budyko.vector, BasinObs$Qmm[Ind.Run])))
  of.kgenp <- (1-(hydroGOF::KGEnp(Qo.budyko.vector, BasinObs$Qmm[Ind.Run])))
  of.lnse <- (abs(hydroGOF::NSE(log(Qo.budyko.vector), log(BasinObs$Qmm[Ind.Run]))))
  of.test01 <- 0.05*(of.nse) + 0.95*(of.nse)
  
  # sumQ <- sum(budyko.out)
  # df.sum <- data.frame(c(sum(budyko.out)/sum(Qobs.budyko)))
  # names(df.sum) <- c("budyko_out")
  # return(df.sum)
  
  # Function return is defined
  return(of.kge)
}

# An initial simulation is executed
1- f.budyko(parm = parm.bud.INI)

# Optimization by GenSA
CNT.budyko <- GenSA(par = parm.bud.INI,
                    fn = f.budyko, 
                    lower = parm.bud.LOW,
                    upper = parm.bud.HIG,
                    control = list(maxit = 1000,
                                   nb.stop.improvement = 500,
                                   max.call = 20000,
                                   verbose = TRUE))

# Optimum budyko parameters are requested
CNT.budyko$par
CNT.budyko$counts

# Final optimum parameters data.frame is created
df.budyko.final.par <- data.frame(as.list(CNT.budyko$par), as.list(CNT.budyko$counts)) 

# Budyko input objects are created
P.budyko <- c(BasinObs$P[Ind.Run])
PE.budyko <- c(BasinObs$E[Ind.Run])
prueba.P <- data.frame(as.list(P.budyko))
prueba.PET <- data.frame(as.list(PE.budyko))
prueba.P <- unname(prueba.P)
prueba.PET <- unname(prueba.PET)

# Budyko model is rerun using optimal parameters
budyko.out <- f.budyko.local(prueba.P, 
                             prueba.PET,
                             g_v = CNT.budyko$par[1],
                             s_v = CNT.budyko$par[2],
                             alpha1_v = CNT.budyko$par[3],
                             alpha2_v = CNT.budyko$par[4],
                             smax_v = CNT.budyko$par[5],
                             d_v = CNT.budyko$par[6])

# Model outputs are isolated
v.budyko.out <- as.vector(budyko.out$q_total)

# Final KGE/NSE are requested
hydroGOF::KGE(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::NSE(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::R2(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::pbias(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::rPearson(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::d(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::md(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal])

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((v.budyko.out - BasinObs$Qmm[Ind.Run.Cal])^2, na.rm = TRUE)
N <- length(v.budyko.out)
K <- 6
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.budyko.gensa.cal <- as.data.frame(t(c(hydroGOF::KGE(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              hydroGOF::NSE(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              hydroGOF::pbias(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              hydroGOF::rPearson(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              hydroGOF::d(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              hydroGOF::md(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              AIC.i,
                                              BIC.i,
                                              hydroGOF::rsr(v.budyko.out, BasinObs$Qmm[Ind.Run.Cal]),
                                              c("budyko"),
                                              c("gensa"),
                                              c("cal"))))

# data.frame variables are renamed
names(df.cont.budyko.gensa.cal) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic","rsr", "model", "method", "stage")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# Calibration run period is defined
Ind.Run = Ind.Run.Val

# Budyko input objects are created
P.budyko <- c(BasinObs$P[Ind.Run])
PE.budyko <- c(BasinObs$E[Ind.Run])
prueba.P <- data.frame(as.list(P.budyko))
prueba.PET <- data.frame(as.list(PE.budyko))
prueba.P <- unname(prueba.P)
prueba.PET <- unname(prueba.PET)

# Budyko model is rerun using optimal parameters
budyko.val.out <- f.budyko.local(prueba.P, 
                                 prueba.PET,
                                 g_v = CNT.budyko$par[1],
                                 s_v = CNT.budyko$par[2],
                                 alpha1_v = CNT.budyko$par[3],
                                 alpha2_v = CNT.budyko$par[4],
                                 smax_v = CNT.budyko$par[5],
                                 d_v = CNT.budyko$par[6])

# Model outputs are isolated
v.budyko.val.out <- as.vector(budyko.val.out$q_total)

# Final KGE/NSE are requested
hydroGOF::KGE(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::NSE(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::R2(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::pbias(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::rPearson(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::d(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])
hydroGOF::md(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val])

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((v.budyko.val.out - BasinObs$Qmm[Ind.Run.Val])^2, na.rm = TRUE)
N <- length(v.budyko.val.out)
K <- 6
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.budyko.gensa.val <- as.data.frame(t(c(hydroGOF::KGE(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              hydroGOF::NSE(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              hydroGOF::pbias(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              hydroGOF::rPearson(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              hydroGOF::d(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              hydroGOF::md(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              AIC.i,
                                              BIC.i,
                                              hydroGOF::rsr(v.budyko.val.out, BasinObs$Qmm[Ind.Run.Val]),
                                              c("budyko"),
                                              c("gensa"),
                                              c("val"))))

# data.frame variables are renamed
names(df.cont.budyko.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Optimization results data.frame is created
df.budyko.final.ofs.def <- rbind(df.cont.budyko.gensa.cal, df.cont.budyko.gensa.val)

# data.frame variables are converted to numeric
df.budyko.final.ofs.def[1:9] <- lapply(df.budyko.final.ofs.def[1:9], as.numeric)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Australian Water Balance Model (AWBM)
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# SOURCES:
# https://ewater.atlassian.net/wiki/spaces/SD550/pages/521929751/Australian+Water+Balance+Model+AWBM+-+SRG

# View default parameter ranges:
str(hydromad.options("awbm"))

# View(HydroTestData)
data(HydroTestData)

# Catchment data.frame is cloned
df.ini <- BasinObs[,c(2,3,1)]

# data.frame variables are renamed
names(df.ini) <- c("P","E","Q")

# New daily dates are assigned to object
inds.cal <- seq(as.Date("2000-01-01 00:00:00"), as.Date("2000-12-25 00:00:00"), by = "day")

# data.frame is transformed to time-series objects
time_series_xts_cal <- xts(df.ini, order.by=as.Date (inds.cal))
time_series_zoo_cal <- zoo(df.ini, order.by=as.Date (inds.cal))

# A AWBM test object is created
test.awbm <- awbm.sim(DATA = time_series_zoo_cal[Ind.Run.Cal],
                      cap.ave = 605,
                      etmult = 435,
                      cap1 = 45,
                      cap2 = 460,
                      cap3 = 923)

# A AWBM model object is created
modx.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "awbm",          # "cwi", "awbm"
                     routing = "expuh",     # "expuh" "powuh"
                     #area1 = 0.134,        # fixed
                     #area2 = 0.433,        # fixed
                     #area3 = 0.433,        # fixed
                     cap.ave = c(40, 2000), # solves C1, C2 and C3
                     etmult = c(0.001, 10),  # Et0 parameter
                     S1_0 = c(0.001, 500),   # quasi-parameter
                     S2_0 = c(0.001, 500),   # quasi-parameter
                     S3_0 = c(0.001, 500),   # quasi-parameter
                     v_s = c(0.001, 1),      # BFI
                     tau_s = c(0.001, 1),    # Kb
                     tau_q = c(0.001, 2))    # Ks

# AWBM model initial state is requested
modx.cal

# AWBM model is optimized
fitx.cal <- fitByOptim(modx.cal,
                       method = "L-BFGS-B",
                       objective = hmadstat("KGE"),
                       samples = 500)

# AWBM optimum calibrated parameters are requested
fitx.cal

# AWBM model performance is requested
summary(fitx.cal)

# a simple grpah is reated
xyplot(fitx.cal, with.P = TRUE, type = c("l", "g"))

# AWBM model is rerun with optimum parameters
mod1.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "awbm", # "cwi", "awbm"
                     routing = "expuh", #"expuh" "powuh"
                     #area1 = 0.134, fixed
                     #area2 = 0.433, fixed
                     #area3 = 0.433, fixed
                     cap.ave = fitx.cal$parlist$cap.ave,
                     etmult = fitx.cal$parlist$etmult,
                     S1_0 = fitx.cal$parlist$S1_0,
                     S2_0 = fitx.cal$parlist$S2_0,
                     S3_0 = fitx.cal$parlist$S3_0,
                     v_s = fitx.cal$parlist$v_s,
                     tau_s = fitx.cal$parlist$tau_s,
                     tau_q = fitx.cal$parlist$tau_q)

# AWBM model prediction are requested
testQ.cal <- predict(mod1.cal, return_state = TRUE)

# A observation vs similations data.frame is created
df.awbm.cal <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], awbm = testQ.cal))

# # A data.frame is requested
# View(df.awbm.cal)

# Performance metrics are requested
hydroGOF::NSE(df.awbm.cal$awbm.X, df.awbm.cal$Q)
hydroGOF::KGE(df.awbm.cal$awbm.X, df.awbm.cal$Q)
hydroGOF::R2(df.awbm.cal$awbm.X, df.awbm.cal$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.awbm.cal$awbm.X - df.awbm.cal$Q)^2, na.rm = TRUE)
N <- length(df.awbm.cal$awbm.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# AWBM final calibrated parameters re requested
area2 = 0.433   # fix
area3 = 0.433   # fix
area1 = 1 - area2 - area3 # fix
cap1 = 0.01 * fitx.cal$parlist$cap.ave/area1
cap2 = 0.33 * fitx.cal$parlist$cap.ave/area2
cap3 = 0.66 * fitx.cal$parlist$cap.ave/area3
tau_s = fitx.cal$parlist$tau_s
v_s = fitx.cal$parlist$v_s
tau_q = fitx.cal$parlist$tau_q
etmult = fitx.cal$parlist$etmult
cap.ave = fitx.cal$parlist$cap.ave

#--------------------------------------------------------------------------------------------------------------
# AWBM hydrological model function is defined for optimization

# Initial fixed parameters are defined
parMM.INI <- c(cap.ave =  100,
               etmult = 1,
               S1_0 = 250, 
               S2_0 = 250, 
               S3_0 = 250, 
               v_s = 0.50,
               tau_s = 0.50,
               tau_q = 0.50)

# Parameter ranges are defined
parMM.LOW <- c(cap.ave =  40,
               etmult = 0.001,
               S1_0 = 0.001, 
               S2_0 = 0.001, 
               S3_0 = 0.001, 
               v_s = 0.001,
               tau_s = 0.001,
               tau_q = 0.001)

parMM.HIG <- c(cap.ave =  2000,
               etmult = 10,
               S1_0 = 500,
               S2_0 = 500,
               S3_0 = 500,
               v_s = 1,
               tau_s = 1,
               tau_q = 2)

# # Sensitivity parameter ranges are defined
# parRanges <- data.frame(min = c(0.1, 260, 0.1, 0.1, 10, 10),
#                         max = c(1.0, 1900, 1.0, 1.0, 100, 60))
# 
# # data.frames variales are defined
# rownames(parRanges) <- c("a", "b", "c", "d", "S_ini", "G_ini")

#--------------------------------------------------------------------------------------------------------------
# AWBM hydrological model function is defined for optimization
AWBM.opti <- function(Param) {
  
  # AWBM model is rerun with optimum parameters
  mod1.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                       warmup = 0,
                       sma = "awbm",        # "cwi", "awbm"
                       routing = "expuh",   # "expuh" "powuh"
                       cap.ave = Param[1],  # solves C1, C2 and C3
                       etmult = Param[2],   # Et0 parameter
                       S1_0 = Param[3],     # quasi-parameter
                       S2_0 = Param[4],     # quasi-parameter
                       S3_0 = Param[5],     # quasi-parameter
                       v_s = Param[6],      # BFI
                       tau_s = Param[7],    # Kb
                       tau_q = Param[8])    # Ks
  
  # AWBM model prediction are requested
  testQ.cal <- predict(mod1.cal, return_state = TRUE)
  
  # A observation vs similations data.frame is created
  df.awbm.cal <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], awbm = testQ.cal))
  
  # A data.frame is requested
  #View(df.awbm.cal)
  
  # Performance metrics are requested
  hydroGOF::NSE(df.awbm.cal$awbm.X, df.awbm.cal$Q)
  hydroGOF::KGE(df.awbm.cal$awbm.X, df.awbm.cal$Q)
  hydroGOF::R2(df.awbm.cal$awbm.X, df.awbm.cal$Q)
  
  # Function return is defined
  return(as.vector(1- hydroGOF::KGE(df.awbm.cal$awbm.X, df.awbm.cal$Q)))
}
#--------------------------------------------------------------------------------------------------------------

# Initial test
AWBM.opti(parMM.INI)

# Optimization by GenSA
CNT.awbm <- GenSA(par = parMM.INI,
                  fn = AWBM.opti, 
                  lower = parMM.LOW,
                  upper = parMM.HIG,
                  control = list(maxit = 1000,
                                 nb.stop.improvement = 500,
                                 max.call = 20000,
                                 verbose = TRUE))

# Optimum awbm parameters are requested
CNT.awbm$par
CNT.awbm$counts

# Final optimum parameters data.frame is created
df.awbm.final.par <- data.frame(as.list(CNT.awbm$par), as.list(CNT.awbm$counts))

# AWBM model is rerun with optimum parameters
mod2.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "awbm", # "cwi", "awbm"
                     routing = "expuh", #"expuh" "powuh"
                     cap.ave = CNT.awbm$par[1],
                     etmult = CNT.awbm$par[2],
                     S1_0 = CNT.awbm$par[3],
                     S2_0 = CNT.awbm$par[4],
                     S3_0 = CNT.awbm$par[5],
                     v_s = CNT.awbm$par[6],
                     tau_s = CNT.awbm$par[7],
                     tau_q = CNT.awbm$par[8])

# AWBM model prediction are requested
testQ.cal2 <- predict(mod2.cal, return_state = TRUE)

# A observation vs similations data.frame is created
df.awbm.cal2 <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], awbm = testQ.cal2))

# A data.frame is requested
# View(df.awbm.cal2)

# Performance metrics are requested
hydroGOF::KGE(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)
hydroGOF::NSE(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)
hydroGOF::R2(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)
hydroGOF::pbias(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)
hydroGOF::rPearson(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)
hydroGOF::d(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)
hydroGOF::md(df.awbm.cal2$awbm.X, df.awbm.cal2$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.awbm.cal2$awbm.X - df.awbm.cal2$Q)^2, na.rm = TRUE)
N <- length(df.awbm.cal2$awbm.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.awbm.gensa.cal <- as.data.frame(t(c(hydroGOF::KGE(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            hydroGOF::NSE(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            hydroGOF::pbias(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            hydroGOF::rPearson(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            hydroGOF::d(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            hydroGOF::md(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            AIC.i,
                                            BIC.i,
                                            hydroGOF::rsr(df.awbm.cal2$awbm.X, df.awbm.cal2$Q),
                                            c("awbm"),
                                            c("gensa"),
                                            c("cal"))))

# data.frame variables are renamed
names(df.cont.awbm.gensa.cal) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# AWBM model is rerun with optimum parameters
mod2.val <- hydromad(time_series_zoo_cal[Ind.Run.Val],
                     warmup = 0,
                     sma = "awbm", # "cwi", "awbm"
                     routing = "expuh", #"expuh" "powuh"
                     cap.ave = CNT.awbm$par[1],
                     etmult = CNT.awbm$par[2],
                     S1_0 = CNT.awbm$par[3],
                     S2_0 = CNT.awbm$par[4],
                     S3_0 = CNT.awbm$par[5],
                     v_s = CNT.awbm$par[6],
                     tau_s = CNT.awbm$par[7],
                     tau_q = CNT.awbm$par[8])

# AWBM model prediction are requested
testQ.val2 <- predict(mod2.val, return_state = TRUE)

# A observation vs similations data.frame is created
df.awbm.val2 <- as.data.frame(cbind(df.ini[Ind.Run.Val,], awbm = testQ.val2))

# A data.frame is requested
#View(df.awbm.val2)

# Performance metrics are requested
hydroGOF::KGE(df.awbm.val2$awbm.X, df.awbm.val2$Q)
hydroGOF::NSE(df.awbm.val2$awbm.X, df.awbm.val2$Q)
hydroGOF::R2(df.awbm.val2$awbm.X, df.awbm.val2$Q)
hydroGOF::pbias(df.awbm.val2$awbm.X, df.awbm.val2$Q)
hydroGOF::rPearson(df.awbm.val2$awbm.X, df.awbm.val2$Q)
hydroGOF::d(df.awbm.val2$awbm.X, df.awbm.val2$Q)
hydroGOF::md(df.awbm.val2$awbm.X, df.awbm.val2$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.awbm.val2$awbm.X - df.awbm.val2$Q)^2, na.rm = TRUE)
N <- length(df.awbm.val2$awbm.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.awbm.gensa.val <- as.data.frame(t(c(hydroGOF::KGE(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            hydroGOF::NSE(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            hydroGOF::pbias(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            hydroGOF::rPearson(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            hydroGOF::d(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            hydroGOF::md(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            AIC.i,
                                            BIC.i,
                                            hydroGOF::rsr(df.awbm.val2$awbm.X, df.awbm.val2$Q),
                                            c("awbm"),
                                            c("gensa"),
                                            c("val"))))

# data.frame variables are renamed
names(df.cont.awbm.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Optimization results data.frame is created
df.awbm.final.ofs.def <- rbind(df.cont.awbm.gensa.cal, df.cont.awbm.gensa.val)

# data.frame variables are converted to numeric
df.awbm.final.ofs.def[1:9] <- lapply(df.awbm.final.ofs.def[1:9], as.numeric)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: IHACRES Catchment Wetness Index (CWI) model
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# SOURCES:
# https://ewater.atlassian.net/wiki/spaces/SD550/pages/521929751/Australian+Water+Balance+Model+AWBM+-+SRG

# A IHACRES.CWI model object is created
modx.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "cwi", # "cwi", "awbm"
                     routing = "expuh", #"expuh" "powuh"
                     #t_ref = c(20, 35),
                     tw = c(0.001, 100),
                     f = c(0.001, 8),
                     scale = c(0.001, 1),
                     s_0 = c(0.001, 1000),
                     v_s = c(0.001, 1),   # Fractional volume for slow flow
                     tau_s = c(0.001, 1), # Time constant for slow flow store [days]
                     tau_q = c(0.001, 2)) # Time constant for quick flow store [days]

# IHACRES.CWI model initial state is requested
modx.cal

# IHACRES.CWI model is optimized
fitx.cal <- fitByOptim(modx.cal,
                       method = "L-BFGS-B",
                       objective = hmadstat("KGE"),
                       samples = 500)

# IHACRES.CWI optimum calibrated parameters are requested
fitx.cal

# IHACRES.CWI model performance is requested
summary(fitx.cal)

# a simple grpah is reated
xyplot(fitx.cal, with.P = TRUE, type = c("l", "g"))

# IHACRES.CWI model is rerun with optimum parameters
mod1.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "cwi", # "cwi", "awbm"
                     routing = "expuh", #"expuh" "powuh"
                      #t_ref = fitx.cal$parlist$t_ref,
                     tw = fitx.cal$parlist$tw,
                     f = fitx.cal$parlist$f,
                     scale = fitx.cal$parlist$scale,
                     s_0 = fitx.cal$parlist$s_0,
                     v_s = fitx.cal$parlist$v_s,
                     tau_s = fitx.cal$parlist$tau_s,
                     tau_q = fitx.cal$parlist$tau_q)

# IHACRES.CWI model prediction are requested
testQ.cal <- predict(mod1.cal, return_state = TRUE)

# A observation vs similations data.frame is created
df.ihacres.cal <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], IHACRES = testQ.cal))

# A data.frame is requested
# View(df.ihacres.cal)

# Performance metrics are requested
hydroGOF::NSE(df.ihacres.cal$IHACRES.X, df.ihacres.cal$Q)
hydroGOF::KGE(df.ihacres.cal$IHACRES.X, df.ihacres.cal$Q)
hydroGOF::R2(df.ihacres.cal$IHACRES.X, df.ihacres.cal$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.ihacres.cal$IHACRES.X - df.ihacres.cal$Q)^2, na.rm = TRUE)
N <- length(df.ihacres.cal$IHACRES.X)
K <- 7
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# IHACRES.CWI final calibrated parameters re requested
tw = fitx.cal$parlist$tw
f = fitx.cal$parlist$f
scale = fitx.cal$parlist$scale
l = fitx.cal$parlist$l
p = fitx.cal$parlist$p
t_ref = fitx.cal$parlist$t_ref
tau_s = fitx.cal$parlist$tau_s
v_s = fitx.cal$parlist$v_s
tau_q = fitx.cal$parlist$tau_q

#--------------------------------------------------------------------------------------------------------------
# IHACRES.CWI hydrological model function is defined for optimization

# Initial fixed parameters are defined
parHH.INI <- c(25,
               5,
               0.0001,
               750, 
               0.50,
               0.50,
               0.50)

# Parameter ranges are defined
parHH.LOW <- c(0.0001,
               0.0001,
               0.0001,
               0.0001,
               0.0001,
               0.0001,
               0.0001)

parHH.HIG <- c(100,
               8,
               1,
               1000,
               1,
               1,
               2)

#--------------------------------------------------------------------------------------------------------------
# IHACRES hydrological model function is defined for optimization
IHACRES.opti <- function(ParamH) {

  # IHACRES model is rerun with optimum parameters
  mod1.cal.hh <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                          warmup = 0,
                          sma = "cwi", # "cwi", "awbm"
                          routing = "expuh", #"expuh" "powuh"
                          #t_ref = c(20, 35),
                          tw = ParamH[1],
                          f = ParamH[2],
                          scale = ParamH[3],
                          s_0 = ParamH[4],
                          v_s = ParamH[5],
                          tau_s = ParamH[6],
                          tau_q = ParamH[7])
  
   # IHACRES model prediction are requested
  testQ.cal <- predict(mod1.cal.hh, return_state = TRUE)
  
  # A observation vs similations data.frame is created
  df.ihacres.cal <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], ihacres = testQ.cal))
  
  # A data.frame is requested
  #View(df.ihacres.cal)
  
  # Performance metrics are requested
  hydroGOF::NSE(df.ihacres.cal$ihacres.X, df.ihacres.cal$Q)
  hydroGOF::KGE(df.ihacres.cal$ihacres.X, df.ihacres.cal$Q)
  hydroGOF::R2(df.ihacres.cal$ihacres.X, df.ihacres.cal$Q)
  
  # Function return is defined
  return(as.vector(1- hydroGOF::KGE(df.ihacres.cal$ihacres.X, df.ihacres.cal$Q)))
}
#--------------------------------------------------------------------------------------------------------------

# Initial test
IHACRES.opti(parHH.INI)

# Optimization by GenSA
CNT.ihacres <- GenSA(par = parHH.INI,
                     fn = IHACRES.opti, 
                     lower = parHH.LOW,
                     upper = parHH.HIG,
                     control = list(maxit = 1000,
                                    nb.stop.improvement = 500,
                                    max.call = 20000,
                                    verbose = TRUE))

# Optimum ihacres parameters are requested
CNT.ihacres$par
CNT.ihacres$counts

# Final optimum parameters data.frame is created
df.ihacres.final.par <- data.frame(as.list(CNT.ihacres$par), as.list(CNT.ihacres$counts)) 

# data.frame variables are renamed
names(df.ihacres.final.par) <- c("tw", "f", "scale", "s_0", "v_s", "tau_s", "tau_q", "calls")

# IHACRES model is rerun with optimum parameters
mod2.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "cwi", # "cwi", "ihacres"
                     routing = "expuh", #"expuh" "powuh"
                     tw = CNT.ihacres$par[1],
                     f = CNT.ihacres$par[2],
                     scale = CNT.ihacres$par[3],
                     s_0 = CNT.ihacres$par[4],
                     v_s = CNT.ihacres$par[5],
                     tau_s = CNT.ihacres$par[6],
                     tau_q = CNT.ihacres$par[7])

# IHACRES model prediction are requested
testQ.cal2 <- predict(mod2.cal, return_state = TRUE)

# A observation vs similations data.frame is created
df.ihacres.cal2 <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], ihacres = testQ.cal2))

# Performance metrics are requested
hydroGOF::KGE(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)
hydroGOF::NSE(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)
hydroGOF::R2(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)
hydroGOF::pbias(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)
hydroGOF::rPearson(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)
hydroGOF::d(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)
hydroGOF::md(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.ihacres.cal2$ihacres.X - df.ihacres.cal2$Q)^2, na.rm = TRUE)
N <- length(df.ihacres.cal2$ihacres.X)
K <- 7
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.ihacres.gensa.cal <- as.data.frame(t(c(hydroGOF::KGE(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               hydroGOF::NSE(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               hydroGOF::pbias(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               hydroGOF::rPearson(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               hydroGOF::d(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               hydroGOF::md(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               AIC.i,
                                               BIC.i,
                                               hydroGOF::rsr(df.ihacres.cal2$ihacres.X, df.ihacres.cal2$Q),
                                               c("ihacres"),
                                               c("gensa"),
                                               c("cal"))))

# data.frame variables are renamed
names(df.cont.ihacres.gensa.cal) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# IHACRES model is rerun with optimum parameters
mod3.val <- hydromad(time_series_zoo_cal[Ind.Run.Val],
                     warmup = 0,
                     sma = "cwi", # "cwi", "ihacres"
                     routing = "expuh", #"expuh" "powuh"
                     tw = CNT.ihacres$par[1],
                     f = CNT.ihacres$par[2],
                     svale = CNT.ihacres$par[3],
                     s_0 = CNT.ihacres$par[4],
                     v_s = CNT.ihacres$par[5],
                     tau_s = CNT.ihacres$par[6],
                     tau_q = CNT.ihacres$par[7])

# IHACRES model prediction are requested
testQ.val2 <- predict(mod3.val, return_state = TRUE)

# A observation vs similations data.frame is created
df.ihacres.val2 <- as.data.frame(cbind(df.ini[Ind.Run.Val,], ihacres = testQ.val2))

# Performance metrics are requested
hydroGOF::KGE(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)
hydroGOF::NSE(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)
hydroGOF::R2(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)
hydroGOF::pbias(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)
hydroGOF::rPearson(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)
hydroGOF::d(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)
hydroGOF::md(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.ihacres.val2$ihacres.X - df.ihacres.val2$Q)^2, na.rm = TRUE)
N <- length(df.ihacres.val2$ihacres.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.ihacres.gensa.val <- as.data.frame(t(c(hydroGOF::KGE(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               hydroGOF::NSE(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               hydroGOF::pbias(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               hydroGOF::rPearson(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               hydroGOF::d(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               hydroGOF::md(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               AIC.i,
                                               BIC.i,
                                               hydroGOF::rsr(df.ihacres.val2$ihacres.X, df.ihacres.val2$Q),
                                               c("ihacres"),
                                               c("gensa"),
                                               c("val"))))

# data.frame variables are renamed
names(df.cont.ihacres.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Optimization results data.frame is created
df.ihacres.final.ofs.def <- rbind(df.cont.ihacres.gensa.cal, df.cont.ihacres.gensa.val)

# data.frame variables are converted to numeric
df.ihacres.final.ofs.def[1:9] <- lapply(df.ihacres.final.ofs.def[1:9], as.numeric)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Simple-Hydrology Model (SimHyd)
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A simhyd model object is created
modv.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                     warmup = 0,
                     sma = "simhyd", # "cwi", "simhyd"
                     routing = NULL, #"expuh" "powuh"
                     INSC = c(0.5, 5),
                     COEFF = c(50, 400),
                     SQ = c(0.0001, 6), #fixed
                     SMSC = c(50, 500), #fixed
                     SUB = c(0.0001, 1), #fixed
                     CRAK = c(0.0001, 1),   # BFI
                     K = c(0.0001, 1), # Ks
                     etmult = c(0.01, 1),
                     GWt0 = c(0.0001, 100),
                     SMSt0 = c(0.0001, 100))
                     #v_s = c(0.001, 1),   # Fractional volume for slow flow
                     #tau_s = c(0.001, 1), # Time constant for slow flow store [days]
                     #tau_q = c(0.001, 2)) # Time constant for quick flow store [days]

# simhyd model initial state is requested
modv.cal

# simhyd model is optimized
fitv.cal <- fitByOptim(modv.cal,
                       method = "L-BFGS-B",
                       objective = hmadstat("KGE"),
                       samples = 1000)

# simhyd optimum calibrated parameters are requested
fitv.cal

# simhyd model performance is requested
summary(fitv.cal)

# a simple grpah is reated
xyplot(fitv.cal, with.P = TRUE, type = c("l", "g"))

# simhyd model is rerun with optimum parameters
modv1.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                      warmup = 0,
                      sma = "simhyd", # "cwi", "simhyd"
                      routing = NULL, #"expuh" "powuh"
                      INSC = fitv.cal$parlist$INSC,
                      COEFF = fitv.cal$parlist$COEFF,
                      SQ = fitv.cal$parlist$SQ, #fixed
                      SMSC = fitv.cal$parlist$SMSC, #fixed
                      SUB = fitv.cal$parlist$SUB, #fixed
                      CRAK = fitv.cal$parlist$CRAK,   # BFI
                      K = fitv.cal$parlist$K,
                      etmult = fitv.cal$parlist$etmult,
                      GWt0 = fitv.cal$parlist$GWt0,
                      SMSt0 = fitv.cal$parlist$SMSt0)
                      #v_s = fitv.cal$parlist$v_s,
                      #tau_s = fitv.cal$parlist$tau_s,
                      #tau_q = fitv.cal$parlist$tau_q) # Ks


# simhyd model prediction are requested
testvQ.cal <- predict(modv1.cal, return_state = TRUE)

# A observation vs similations data.frame is created
df.simhyd.cal <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], simhyd = testvQ.cal))

# # A data.frame is requested
# View(df.simhyd.cal)

# Performance metrics are requested
hydroGOF::NSE(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)
hydroGOF::KGE(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)
hydroGOF::R2(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.simhyd.cal$simhyd.X - df.simhyd.cal$Q)^2, na.rm = TRUE)
N <- length(df.simhyd.cal$simhyd.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# SimHyd hydrological model function is defined for optimization

# Initial fixed parameters are defined
parSS.INI <- c(INSC =  2,
               COEFF = 200,
               SQ = 0.1, 
               SMSC = 280, 
               SUB = 0.5, 
               CRAK = 0.9,
               K = 0.1,
               etmult = 0.1,
               GWt0 = 60,
               SMSt0 = 60)

# Parameter ranges are defined
parSS.LOW <- c(INSC =  0.50,
               COEFF = 50.0,
               SQ = 0.000001, 
               SMSC = 50.0, 
               SUB = 0.0001, 
               CRAK = 0.0001,
               K = 0.0001,
               etmult = 0.0001,
               GWt0 = 0.0001,
               SMSt0 = 0.0001)

parSS.HIG <- c(INSC =  5.0,
               COEFF = 400.0,
               SQ = 6.00,
               SMSC = 500.0,
               SUB = 1.0,
               CRAK = 1.0,
               K = 1.0,
               etmult = 1.0,
               GWt0 = 100.00,
               SMSt0 = 100.00)

#--------------------------------------------------------------------------------------------------------------
# simhyd hydrological model function is defined for optimization
simhyd.opti <- function(Param) {
  
  # simhyd model is rerun with optimum parameters
  mods1.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                        warmup = 0,
                        sma = "simhyd", # "cwi", "simhyd"
                        routing = NULL, #"expuh" "powuh"
                        INSC = Param[1],
                        COEFF = Param[2],
                        SQ = Param[3], #fixed
                        SMSC = Param[4], #fixed
                        SUB = Param[5], #fixed
                        CRAK = Param[6],   # BFI
                        K = Param[7],
                        etmult = Param[8],
                        GWt0 = Param[9],
                        SMSt0 = Param[10])
                        #v_s = fitv.cal$parlist$v_s,
                        #tau_s = fitv.cal$parlist$tau_s,
                        #tau_q = fitv.cal$parlist$tau_q) # Ks
  
  # simhyd model prediction are requested
  testQS.cal <- predict(mods1.cal, return_state = TRUE)
  
  # A observation vs similations data.frame is created
  df.simhyd.cal <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], simhyd = testQS.cal))
  
  # A data.frame is requested
  #View(df.simhyd.cal)
  
  # Performance metrics are requested
  hydroGOF::NSE(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)
  hydroGOF::KGE(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)
  hydroGOF::R2(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)
  
  # Function return is defined
  return(as.vector(1- hydroGOF::KGE(df.simhyd.cal$simhyd.U, df.simhyd.cal$Q)))
}
#--------------------------------------------------------------------------------------------------------------

# Initial test
simhyd.opti(parSS.INI)

# Optimization by GenSA
CNT.simhyd <- GenSA(par = parSS.INI,
                    fn = simhyd.opti, 
                    lower = parSS.LOW,
                    upper = parSS.HIG,
                    control = list(maxit = 1000,
                                   nb.stop.improvement = 500,
                                   max.call = 20000,
                                   verbose = TRUE))

# Optimum simhyd parameters are requested
CNT.simhyd$par
CNT.simhyd$counts

# Final optimum parameters data.frame is created
df.simhyd.final.par <- data.frame(as.list(CNT.simhyd$par), as.list(CNT.simhyd$counts))

# simhyd model is rerun with optimum parameters
modS2.cal <- hydromad(time_series_zoo_cal[Ind.Run.Cal],
                      warmup = 0,
                      sma = "simhyd", # "cwi", "simhyd"
                      routing = "expuh", #"expuh" "powuh"
                      INSC = CNT.simhyd$par[1],
                      COEFF = CNT.simhyd$par[2],
                      SQ = CNT.simhyd$par[3],
                      SMSC = CNT.simhyd$par[4],
                      SUB = CNT.simhyd$par[5],
                      CRAK = CNT.simhyd$par[6],
                      K = CNT.simhyd$par[7],
                      etmult = CNT.simhyd$par[8],
                      GWt0 = CNT.simhyd$par[9],
                      SMSt0 = CNT.simhyd$par[10])

# simhyd model prediction are requested
testQS.cal2 <- predict(modS2.cal, return_state = TRUE)

# A observation vs similations data.frame is created
df.simhyd.cal2 <- as.data.frame(cbind(df.ini[Ind.Run.Cal,], simhyd = testQS.cal2))

# A data.frame is requested
# View(df.simhyd.cal2)

# Performance metrics are requested
hydroGOF::KGE(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)
hydroGOF::NSE(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)
hydroGOF::R2(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)
hydroGOF::pbias(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)
hydroGOF::rPearson(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)
hydroGOF::d(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)
hydroGOF::md(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.simhyd.cal2$simhyd.X - df.simhyd.cal2$Q)^2, na.rm = TRUE)
N <- length(df.simhyd.cal2$simhyd.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------


# A data.frame container is defined
df.cont.simhyd.gensa.cal <- as.data.frame(t(c(hydroGOF::KGE(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              hydroGOF::NSE(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              hydroGOF::pbias(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              hydroGOF::rPearson(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              hydroGOF::d(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              hydroGOF::md(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              AIC.i,
                                              BIC.i,
                                              hydroGOF::rsr(df.simhyd.cal2$simhyd.X, df.simhyd.cal2$Q),
                                              c("simhyd"),
                                              c("gensa"),
                                              c("cal"))))

# data.frame variables are renamed
names(df.cont.simhyd.gensa.cal) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# simhyd model is rerun with optimum parameters
modS2.val <- hydromad(time_series_zoo_cal[Ind.Run.Val],
                      warmup = 0,
                      sma = "simhyd", # "cwi", "simhyd"
                      routing = "expuh", #"expuh" "powuh"
                      INSC = CNT.simhyd$par[1],
                      COEFF = CNT.simhyd$par[2],
                      SQ = CNT.simhyd$par[3],
                      SMSC = CNT.simhyd$par[4],
                      SUB = CNT.simhyd$par[5],
                      CRAK = CNT.simhyd$par[6],
                      K = CNT.simhyd$par[7],
                      etmult = CNT.simhyd$par[8],
                      GWt0 = CNT.simhyd$par[9],
                      SMSt0 = CNT.simhyd$par[10])

# simhyd model prediction are requested
testSQ.val2 <- predict(modS2.val, return_state = TRUE)

# A observation vs similations data.frame is created
df.simhyd.val2 <- as.data.frame(cbind(df.ini[Ind.Run.Val,], simhyd = testSQ.val2))

# A data.frame is requested
#View(df.simhyd.val2)

# Performance metrics are requested
hydroGOF::KGE(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)
hydroGOF::NSE(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)
hydroGOF::R2(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)
hydroGOF::pbias(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)
hydroGOF::rPearson(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)
hydroGOF::d(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)
hydroGOF::md(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q)

#------------------------------------------------------------------
# https://itfeature.com/regression/selection/akaike-information-criteria/
# https://itfeature.com
SSe <- sum((df.simhyd.val2$simhyd.X - df.simhyd.val2$Q)^2, na.rm = TRUE)
N <- length(df.simhyd.val2$simhyd.X)
K <- 8
AIC.i <- N*(log(SSe/N)) + 2*(K+1)
print(AIC.i)
BIC.i <- N*(log(SSe/N)) + (log(N))*(K+1)
print(BIC.i)
#------------------------------------------------------------------

# A data.frame container is defined
df.cont.simhyd.gensa.val <- as.data.frame(t(c(hydroGOF::KGE(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              hydroGOF::NSE(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              hydroGOF::pbias(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              hydroGOF::rPearson(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              hydroGOF::d(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              hydroGOF::md(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              AIC.i,
                                              BIC.i,
                                              hydroGOF::rsr(df.simhyd.val2$simhyd.X, df.simhyd.val2$Q),
                                              c("simhyd"),
                                              c("gensa"),
                                              c("val"))))

# data.frame variables are renamed
names(df.cont.simhyd.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Optimization results data.frame is created
df.simhyd.final.ofs.def <- rbind(df.cont.simhyd.gensa.cal, df.cont.simhyd.gensa.val)

# data.frame variables are converted to numeric
df.simhyd.final.ofs.def[1:9] <- lapply(df.simhyd.final.ofs.def[1:9], as.numeric)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Model Ensemble
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#------------------------------------------------------------------------------
# Calibration
#------------------------------------------------------------------------------

# An ensemble data.frame is created
df.ensemble.out <- data.frame(BasinObs$Qmm[Ind.Run.Cal], # observed Q
                              BasinObs$P[Ind.Run.Cal],   # observed P
                              BasinObs$E[Ind.Run.Cal],   # observed ET0
                              df.GR2M.gen$QSIM,          # Q gr2m
                              acbd.out,                  # Q abcd
                              v.budyko.out,              # Q budyko
                              df.awbm.cal2$awbm.X,       # Q awbm
                              df.ihacres.cal2$ihacres.X, # Q ihacres
                              df.simhyd.cal2$simhyd.X)   # Q simhyd

# data.frame variables are renamed
names(df.ensemble.out) <-c("observed_Q",
                           "observed_P",
                           "observed_ET0",
                           "Q_gr2m",
                           "Q_abcd",
                           "Q budyko",
                           "Q_awbm",
                           "Q_ihacres",
                           "Q simhyd")

# The mean value is calculated
df.ensemble.out$mean <- rowMeans(df.ensemble.out[, 4:8])

# Final KGE/NSE are requested
hydroGOF::KGE(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::NSE(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::R2(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::pbias(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::rSpearman(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::d(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::md(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])
hydroGOF::rsr(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal])

# A data.frame container is defined
df.cont.ensemble.gensa.cal <- as.data.frame(t(c(hydroGOF::KGE(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                hydroGOF::NSE(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                hydroGOF::pbias(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                hydroGOF::R2(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                hydroGOF::d(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                hydroGOF::md(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                1000,
                                                1000,
                                                hydroGOF::rsr(df.ensemble.out$mean, BasinObs$Qmm[Ind.Run.Cal]),
                                                c("ensemble"),
                                                c("gensa"),
                                                c("cal"))))

# data.frame variables are renamed
names(df.cont.ensemble.gensa.cal) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Catchment domain is added to data.frame
df.cont.ensemble.gensa.cal$catchment <- c("tempisque")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------

# An ensemble data.frame is created
df.ensemble.out.val <- data.frame(BasinObs$Qmm[Ind.Run.Val], # observed Q
                                  BasinObs$P[Ind.Run.Val],   # observed P
                                  BasinObs$E[Ind.Run.Val],   # observed ET0
                                  df.GR2M.gen.val$QSIM,      # Q gr2m
                                  acbd.out.val,              # Q abcd
                                  v.budyko.val.out,          # Q budyko
                                  df.awbm.val2$awbm.X,       # Q awbm
                                  df.ihacres.val2$ihacres.X, # Q ihacres
                                  df.simhyd.val2$simhyd.X)   # Q simhyd

# data.frame variables are renamed
names(df.ensemble.out.val) <-c("observed_Q",
                               "observed_P",
                               "observed_ET0",
                               "Q_gr2m",
                               "Q_abcd",
                               "Q budyko",
                               "Q_awbm",
                               "Q_ihacres",
                               "Q simhyd")

# The mean value is calculated
df.ensemble.out.val$mean <- rowMeans(df.ensemble.out.val[, 4:8])

# A data.frame container is defined
df.cont.ensemble.gensa.val <- as.data.frame(t(c(hydroGOF::KGE(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                hydroGOF::NSE(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                hydroGOF::pbias(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                hydroGOF::R2(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                hydroGOF::d(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                hydroGOF::md(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                1000,
                                                1000,
                                                hydroGOF::rsr(df.ensemble.out.val$mean, BasinObs$Qmm[Ind.Run.Val]),
                                                c("ensemble"),
                                                c("gensa"),
                                                c("val"))))

# data.frame variables are renamed
names(df.cont.ensemble.gensa.val) <- c("kge", "nse", "pbias", "corr", "d", "md", "aic", "bic", "rsr", "model", "method", "stage")

# Catchment domain is added to data.frame
df.cont.ensemble.gensa.val$catchment <- c("tempisque")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Summary-exploratory data.frames and objects exports
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# GenSA number of model calls are requested
calls.gr2m <- CNT.gr2m$counts
calls.abcd <- CNT.abcd$counts
calls.budyko <- CNT.budyko$counts
calls.awbm <- CNT.awbm$counts
calls.ihacres <- CNT.ihacres$counts
calls.simhyd <- CNT.simhyd$counts

# A data.frame container is defined
df.cont.ensemble.gensa.calls <- as.data.frame(t(c(calls.gr2m,
                                                  calls.abcd,
                                                  calls.budyko,
                                                  calls.awbm,
                                                  calls.ihacres,
                                                  calls.simhyd,
                                                  c("tempisque"))))

# data.frame variables are renamed
names(df.cont.ensemble.gensa.calls) <- c("gr2m",
                                         "abcd",
                                         "budyko",
                                         "awbm",
                                         "ihacres",
                                         "simhyd",
                                         "catchment")

# Catchment domain is added to data.frame
df.gr2m.final.par$catchment <- c("tempisque")
df.abcd.final.par$catchment <- c("tempisque")
df.budyko.final.par$catchment <- c("tempisque")
df.awbm.final.par$catchment <- c("tempisque")
df.ihacres.final.par$catchment <- c("tempisque")
df.simhyd.final.par$catchment <- c("tempisque")

# data.frame is requested
View(df.gr2m.final.par)
View(df.abcd.final.par)
View(df.budyko.final.par)
View(df.awbm.final.par)
View(df.ihacres.final.par)
View(df.simhyd.final.par)

# Catchment domain is added to data.frame
df.gr2m.final.ofs.def$catchment <- c("tempisque")
df.abcd.final.ofs.def$catchment <- c("tempisque")
df.budyko.final.ofs.def$catchment <- c("tempisque")
df.awbm.final.ofs.def$catchment <- c("tempisque")
df.ihacres.final.ofs.def$catchment <- c("tempisque")
df.simhyd.final.ofs.def$catchment <- c("tempisque")

# Optimization results data.frame is created
df.total.final.ofs.def <- rbind(df.gr2m.final.ofs.def,
                                df.abcd.final.ofs.def,
                                df.budyko.final.ofs.def,
                                df.awbm.final.ofs.def,
                                df.ihacres.final.ofs.def,
                                df.simhyd.final.ofs.def,
                                df.cont.ensemble.gensa.cal,
                                df.cont.ensemble.gensa.val)

# Local optimization is deleted
df.total.final.ofs.def <- df.total.final.ofs.def[-c(1), ]

# data.frame is requested
# View(df.gr2m.final.ofs.def)
# View(df.abcd.final.ofs.def)
# View(df.budyko.final.ofs.def)
# View(df.awbm.final.ofs.def)
# View(df.ihacres.final.ofs.def)
# View(df.simhyd.final.ofs.def)
View(df.total.final.ofs.def)

# An output file name is concatenated
t.output01 <- paste("Multi_Hydro_Ensemble_tempisque_SECOND-ROUND",".xlsx", sep = "")

# An export/import XLS object is created
xls.object <- list("df.total.final.ofs.def" = df.total.final.ofs.def, #1
                   "df.gr2m.final.par" = df.gr2m.final.par,           #2
                   "df.abcd.final.par" = df.abcd.final.par,           #3
                   "df.budyko.final.par" = df.budyko.final.par,       #4
                   "df.awbm.final.par" = df.awbm.final.par,           #5
                   "df.ihacres.final.par" = df.ihacres.final.par,     #6
                   "df.simhyd.final.par" = df.simhyd.final.par,       #7
                   "df.ensemble.out.cal" = df.ensemble.out,           #8
                   "df.ensemble.out.val" = df.ensemble.out.val,       #9
                   "df.cont.ensemble.gensa.calls" = df.cont.ensemble.gensa.calls)       #10
write.xlsx(xls.object, file = t.output01, rowNames = TRUE)

# Execution time is calculated
stm <- (proc.time() - ptm)
print((stm[1])/60)

#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------


                                       































# ++++++++++++++++++++++++++++++++++++++++++++++++++++

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Local Sensitivity Analysis
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: GR2M Hydrological model
# FME - Flexible Modelling Environment for Inverse Modelling
#--------------------------------------------------------------------------------------------------------------

# Control parameters are defined
param.control.gr2m <- c(X1 =  1000, X2 = 1.00) # Nominal literature values

# Sensitivity parameter ranges are defined
param.ranges.control <- data.frame(min = c(param.control.gr2m[1]*0.90, param.control.gr2m[2]*0.90),
                                   max = c(param.control.gr2m[1]*1.10, param.control.gr2m[2]*1.10))

# data.frames variales are defined
rownames(param.ranges.control) <- c("X1", "X2")

# GR2M hydrological model is run with control parameters
df.GR2M.out.control <- f.GR2M.prod(param.control.gr2m)

# GR2M hydrological model is run with control parameters
df.GR2M.out.sens <- f.GR2M.sens(param.control.gr2m)

# Sensitivity Analysis via sensRange {FME} for parameter X1
sR.GR2M.X1 <- sensRange(func = f.GR2M.sens,
                        parms = param.control.gr2m, # Param.prod.gensa, # these are the optimum parameters!!!
                        dist = "grid",
                        sensvar = c("sum_obs",     #1
                                    "sum_sim",     #2
                                    "sswr",        #3
                                    "nse",         #4
                                    "lnse",        #5
                                    "kge",         #6
                                    "kgenp",       #7
                                    "pbias",       #8
                                    "ia",          #9
                                    "mia",         #10
                                    "rsr",         #11
                                    "corr",        #12
                                    "rmse",        #13
                                    "of_test_01"), #14
                        parRange = param.ranges.control[c(1),], # parameter X1
                        num = 100)

# Normalized Sensitivity Coefficient (NSC) data.frame is created
df.sR.GR2M.X1.nsc <- sR.GR2M.X1[, c(1,7)]

# sum_obs values are replaced by synthetic control values
df.sR.GR2M.X1.nsc$kge_control <- as.numeric(df.GR2M.out.control$kge)

# X1 optimal value is added to data.frame
df.sR.GR2M.X1.nsc$X1_control <- param.control.gr2m[1]

# NSC SENSAN L2 norm sensitivities are calculated
df.sR.GR2M.X1.nsc$L2 <- (abs(df.sR.GR2M.X1.nsc$kge - df.sR.GR2M.X1.nsc$kge_control)/
                           abs(df.sR.GR2M.X1.nsc$X1 - df.sR.GR2M.X1.nsc$X1_control)) *
  (abs(df.sR.GR2M.X1.nsc$X1_control)/abs(df.sR.GR2M.X1.nsc$kge_control))

# A simple plot is created
plot(df.sR.GR2M.X1.nsc$X1, df.sR.GR2M.X1.nsc$L2)

# A simple boxplot is created
boxplot(df.sR.GR2M.X1.nsc$L2)

# The median NSC value is requested
median(df.sR.GR2M.X1.nsc$L2)

# Sensitivity Analysis via sensRange {FME} for parameter X2
sR.GR2M.X2 <- sensRange(func = f.GR2M.sens,
                        parms = param.control.gr2m, # these are the optimum parameters!!!
                        dist = "grid",
                        sensvar = c("sum_obs",     #1
                                    "sum_sim",     #2
                                    "sswr",        #3
                                    "nse",         #4
                                    "lnse",        #5
                                    "kge",         #6
                                    "kgenp",       #7
                                    "pbias",       #8
                                    "ia",          #9
                                    "mia",         #10
                                    "rsr",         #11
                                    "corr",        #12
                                    "rmse",        #13
                                    "of_test_01"), #14
                        parRange = param.ranges.control[c(2),], # parameter X2
                        num = 100)

# Normalized Sensitivity Coefficient (NSC) data.frame is created
df.sR.GR2M.X2.nsc <- sR.GR2M.X2[, c(1,7)]

# sum_obs values are replaced by synthetic control values
df.sR.GR2M.X2.nsc$kge_control <- as.numeric(df.GR2M.out.control$kge)

# X2 optimal value is added to data.frame
df.sR.GR2M.X2.nsc$X2_control <- param.control.gr2m[2]

# NSC SENSAN L2 norm sensitivities are calculated
df.sR.GR2M.X2.nsc$L2 <- (abs(df.sR.GR2M.X2.nsc$kge - df.sR.GR2M.X2.nsc$kge_control)/
                           abs(df.sR.GR2M.X2.nsc$X2 - df.sR.GR2M.X2.nsc$X2_control)) *
  (abs(df.sR.GR2M.X2.nsc$X2_control)/abs(df.sR.GR2M.X2.nsc$kge_control))

# A simple plot is created
plot(df.sR.GR2M.X2.nsc$X2, df.sR.GR2M.X2.nsc$L2)

# A simple boxplot is created
boxplot(df.sR.GR2M.X2.nsc$L2)

# The median NSC value is requested
median(df.sR.GR2M.X2.nsc$L2)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: ABCD Hydrological model
# FME - Flexible Modelling Environment for Inverse Modelling
#--------------------------------------------------------------------------------------------------------------

# Initial fixed parameters are defined
parm.INI.abcd <- c(a =  0.50, b = 350.00, c = 0.50,  d = 0.50, S_ini = 500.00, G_ini = 60.0)

# Sensitivity parameter ranges are defined
param.ranges.control.abcd <- data.frame(min = c(parm.INI.abcd [1]*0.90,
                                                parm.INI.abcd [2]*0.90,
                                                parm.INI.abcd [3]*0.90,
                                                parm.INI.abcd [4]*0.90,
                                                parm.INI.abcd [5]*0.90,
                                                parm.INI.abcd [6]*0.90),
                                        max = c(parm.INI.abcd [1]*1.10,
                                                parm.INI.abcd [2]*1.10,
                                                parm.INI.abcd [3]*1.10,
                                                parm.INI.abcd [4]*1.10,
                                                parm.INI.abcd [5]*1.10,
                                                parm.INI.abcd [6]*1.10))

# data.frames variales are defined
rownames(param.ranges.control.abcd) <- c("a", "b", "c", "d", "S_ini", "G_ini")

# ABCD hydrological model function is defined for sensitivity analysis
f.ABCD.sens <- function(parm) {
  
  # ABCD input objects are created
  P.abcd <- BasinObs$P[Ind.Run.Cal]
  PE.abcd <- BasinObs$E[Ind.Run.Cal]
  Qobs.abcd <- BasinObs$Qmm[Ind.Run.Cal]
  
  # abcdQest {hydrosystems} function is requested
  acbd.out <- abcdQest(parm = parm[1:4],
                       P = P.abcd,
                       PE = PE.abcd,
                       S_ini = parm[5],
                       G_ini = parm[6],
                       print.all = FALSE)
  
  # ABCD output objects are created
  acbd.out <- as.vector(acbd.out)
  
  # Metric functions are defined
  of.sse <- sum((Qobs.abcd - acbd.out)^2,     na.rm = TRUE)
  of.nse <- (1-(hydroGOF::NSE(acbd.out, Qobs.abcd)))
  of.kge <- (1-(hydroGOF::KGE(acbd.out, Qobs.abcd)))
  of.kgenp <- (1-(hydroGOF::KGEnp(acbd.out, Qobs.abcd)))
  of.lnse <- (abs(hydroGOF::NSE(log(acbd.out), log(Qobs.abcd))))
  of.test01 <- 0.05*(of.nse) + 0.95*(of.nse)
  
  # sumQ <- sum(acbd.out)
  # df.sum <- data.frame(c(sum(acbd.out)/sum(Qobs.abcd)))
  # names(df.sum) <- c("abcd_out")
  # return(df.sum)
  
  # Output data.frame is created
  df.sum <- data.frame(kge = c(of.kge)) 
  
  # Function return is defined
  return(df.sum)
}

# ABCD hydrological model is run with control parameters
df.ABCD.out.control <- f.ABCD.sens(parm.INI.abcd)

# Sensitivity Analysis via sensRange {FME} for parameter A
sR1.ABCD.A <- sensRange(func = f.ABCD.sens,
                 parms = parm.INI.abcd, # these are the control parameters!!!
                 dist = "grid",
                 sensvar = c("kge"),
                 parRange = param.ranges.control.abcd[c(1),], # parameter A
                 num = 100)

# Normalized Sensitivity Coefficient (NSC) data.frame is created
df.sR.ABCD.A.nsc <- sR1.ABCD.A

# sum_obs values are replaced by synthetic control values
df.sR.ABCD.A.nsc$kge_control <- as.numeric(df.ABCD.out.control$kge)

# A optimal value is added to data.frame
df.sR.ABCD.A.nsc$A_control <- parm.INI.abcd[1]

# NSC SENSAN L2 norm sensitivities are calculated
df.sR.ABCD.A.nsc$L2 <- (abs(df.sR.ABCD.A.nsc$kge - df.sR.ABCD.A.nsc$kge_control)/
                          abs(df.sR.ABCD.A.nsc$a - df.sR.ABCD.A.nsc$A_control)) *
  (abs(df.sR.ABCD.A.nsc$A_control)/abs(df.sR.ABCD.A.nsc$kge_control))

# A simple plot is created
plot(df.sR.ABCD.A.nsc$a, df.sR.ABCD.A.nsc$L2)

# A simple boxplot is created
boxplot(df.sR.ABCD.A.nsc$L2)

# The median NSC value is requested
median(df.sR.ABCD.A.nsc$L2)

# Sensitivity Analysis via sensRange {FME} for parameter B
sR1.ABCD.B <- sensRange(func = f.ABCD.sens,
                        parms = parm.INI.abcd, # these are the control parameters!!!
                        dist = "grid",
                        sensvar = c("kge"),
                        parRange = param.ranges.control.abcd[c(2),], # parameter A
                        num = 100)

# Normalized Sensitivity Coefficient (NSC) data.frame is created
df.sR.ABCD.B.nsc <- sR1.ABCD.B

# sum_obs values are replaced by synthetic control values
df.sR.ABCD.B.nsc$kge_control <- as.numeric(df.ABCD.out.control$kge)

# A optimal value is added to data.frame
df.sR.ABCD.B.nsc$B_control <- parm.INI.abcd[2]

# NSC SENSAN L2 norm sensitivities are calculated
df.sR.ABCD.B.nsc$L2 <- (abs(df.sR.ABCD.B.nsc$kge - df.sR.ABCD.B.nsc$kge_control)/
                          abs(df.sR.ABCD.B.nsc$b - df.sR.ABCD.B.nsc$B_control)) *
  (abs(df.sR.ABCD.B.nsc$B_control)/abs(df.sR.ABCD.B.nsc$kge_control))

# A simple plot is created
plot(df.sR.ABCD.B.nsc$b, df.sR.ABCD.B.nsc$L2)

# A simple boxplot is created
boxplot(df.sR.ABCD.B.nsc$L2)

# The median NSC value is requested
median(df.sR.ABCD.B.nsc$L2)

# Sensitivity Analysis via sensRange {FME} for parameter C
sR1.ABCD.C <- sensRange(func = f.ABCD.sens,
                        parms = parm.INI.abcd, # these are the control parameters!!!
                        dist = "grid",
                        sensvar = c("kge"),
                        parRange = param.ranges.control.abcd[c(3),], # parameter A
                        num = 100)

# Normalized Sensitivity Coefficient (NSC) data.frame is created
df.sR.ABCD.C.nsc <- sR1.ABCD.C

# sum_obs values are replaced by synthetic control values
df.sR.ABCD.C.nsc$kge_control <- as.numeric(df.ABCD.out.control$kge)

# A optimal value is added to data.frame
df.sR.ABCD.C.nsc$C_control <- parm.INI.abcd[3]

# NSC SENSAN L2 norm sensitivities are calculated
df.sR.ABCD.C.nsc$L2 <- (abs(df.sR.ABCD.C.nsc$kge - df.sR.ABCD.C.nsc$kge_control)/
                          abs(df.sR.ABCD.C.nsc$c - df.sR.ABCD.C.nsc$C_control)) *
  (abs(df.sR.ABCD.C.nsc$C_control)/abs(df.sR.ABCD.C.nsc$kge_control))

# A simple plot is created
plot(df.sR.ABCD.C.nsc$c, df.sR.ABCD.C.nsc$L2)

# A simple boxplot is created
boxplot(df.sR.ABCD.C.nsc$L2)

# The median NSC value is requested
median(df.sR.ABCD.C.nsc$L2)

# Sensitivity Analysis via sensRange {FME} for parameter C
sR1.ABCD.D <- sensRange(func = f.ABCD.sens,
                        parms = parm.INI.abcd, # these are the control parameters!!!
                        dist = "grid",
                        sensvar = c("kge"),
                        parRange = param.ranges.control.abcd[c(4),], # parameter A
                        num = 100)

# Normalized Sensitivity Coefficient (NSC) data.frame is created
df.sR.ABCD.D.nsc <- sR1.ABCD.D

# sum_obs values are replaced by synthetic control values
df.sR.ABCD.D.nsc$kge_control <- as.numeric(df.ABCD.out.control$kge)

# A optimal value is added to data.frame
df.sR.ABCD.D.nsc$D_control <- parm.INI.abcd[4]

# NSC SENSAN L2 norm sensitivities are calculated
df.sR.ABCD.D.nsc$L2 <- (abs(df.sR.ABCD.D.nsc$kge - df.sR.ABCD.D.nsc$kge_control)/
                          abs(df.sR.ABCD.D.nsc$d - df.sR.ABCD.D.nsc$D_control)) *
  (abs(df.sR.ABCD.D.nsc$D_control)/abs(df.sR.ABCD.D.nsc$kge_control))

# A simple plot is created
plot(df.sR.ABCD.D.nsc$d, df.sR.ABCD.D.nsc$L2)

# A simple boxplot is created
boxplot(df.sR.ABCD.D.nsc$L2)

# The median NSC value is requested
median(df.sR.ABCD.D.nsc$L2)

#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------









# 
# 
# CNT.02 <- DEoptim(fn = f.HBV02, 
#                   lower = Param.LOW,
#                   upper = Param.HIG,
#                   control = DEoptim.control(itermax = 100,
#                                             trace = TRUE))
# 
# CNT.02$optim[c("bestmem","bestval","nfeval", "iter")]




# Grubbs test {outliers} for outliers is performed
#grubbs.X2 <- grubbs.test(df.sR.GR2M.X2.nsc$L2)

# range for 95% QQ is isolated
#qq.X2 <- quantile(df.sR.GR2M.X2.nsc$L2,probs=c(0.975))

# Outliers L2 values are deleted
#df.sR.GR2M.X2.nsc <- subset(df.sR.GR2M.X2.nsc, L2 < qq.X2 )

