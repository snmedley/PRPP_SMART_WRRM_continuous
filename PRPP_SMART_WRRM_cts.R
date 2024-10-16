##### Set-up: Load packages and create function #####
# Load necessary packages
library(LaplacesDemon)
library(geepack)
library(geometry)
library(dplyr)
library(MASS)
library(data.table)
library(broom)
library(msm)
library(kableExtra)

# function to compute confidence interval from geeglm output
confint.geeglm <- function(object, parm, level = 0.95, ...) {
  cc <- coef(summary(object))
  mult <- qnorm((1+level)/2)
  citab <- with(as.data.frame(cc),
                cbind(lwr=Estimate-mult*Std.err,
                      upr=Estimate+mult*Std.err))
  rownames(citab) <- rownames(cc)
  citab[parm,]
}

##### Choose Data Generation Settings #####
# Choose home directory where the data generation file is stored. Files for data generation scenarios
# must be in a folder named "Scenarios" within the home directory.
homedir = "C:/Users/snmed/OneDrive/Documents/GSRA/Project_PCORI_PRPP/ContinuousOutcome/GitHubCode"
# Choose a directory for output files. The directory is set to homedir by default. 
outdir = "C:/Users/snmed/OneDrive/Documents/GSRA/Project_PCORI_PRPP/ContinuousOutcome/GitHubCode/Results"
# Choose sample size
N = 500
# Choose preference/response rate scenario by setting scenario = a, b, or c
scenario = "a" 
setwd(paste0(homedir,"/Scenarios"))
# Load values for pNP_target, pTheta_target, pNP2_target, pTheta2_target, Pa, Pa1, Pb, Pb1
source(paste0("scenario_", scenario, ".R")) 
# Choose preference augmented DTR effect scenario by setting type = 1, 2, 3, or 4 and size = small, moderate, or large
type = 1
size = "small"
# Loads the values for the 20 pathway means
source(paste0("type", type, "_", size, ".R"))
# Load data generation function (depends on chosen model)
setwd(homedir)
source("PRPP_SMART_DataGen_cts.R")
# Choose sigma2 (data variability at trial pathway level) 
sigma2 = 36
# Choose desired number of simulations
n.sim = 500

##### Calculate True Parameter and DTR Values #####
#### DTRs ####
# True DTR values for PRPP-SMART -- depends on Pa, Pa1, Pb, Pb1 from scenario file
# Also depends on pathway means specified in preference setting file
# Does NOT depend on the chosen model
dtr_names = c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
expected_pref <- c()  # expected DTR mean outcomes
expected_pref[1] <- Pa * muRAR + (1 - Pa) * muRANRRC  #AAC00
expected_pref[2] <- Pa * muRAR + (1 - Pa) * muRANRRD  #AAD00
expected_pref[3] <- Pb * muRBR + (1 - Pb) * muRBNRRC  #BBC00
expected_pref[4] <- Pb * muRBR + (1 - Pb) * muRBNRRD  #BBD00
expected_pref[5] <- Pa * muRAR + (1 - Pa) * muRANRPC #AAC01
expected_pref[6] <- Pa * muRAR + (1 - Pa) * muRANRPD #AAD01
expected_pref[7] <- Pb * muRBR + (1 - Pb) * muRBNRPC #BBC01
expected_pref[8] <- Pb * muRBR + (1 - Pb) * muRBNRPD #BBD01
expected_pref[9] <- Pa1 * muPAR + (1 - Pa1) * muPANRRC #AAC10
expected_pref[10] <- Pa1 * muPAR + (1 - Pa1) * muPANRRD #AAD10
expected_pref[11] <- Pb1 * muPBR + (1 - Pb1) * muPBNRRC #BBC10
expected_pref[12] <- Pb1 * muPBR + (1 - Pb1) * muPBNRRD #BBD10
expected_pref[13] <- Pa1 * muPAR + (1 - Pa1) * muPANRPC #AAC11
expected_pref[14] <- Pa1 * muPAR + (1 - Pa1) * muPANRPD #AAD11
expected_pref[15] <- Pb1 * muPBR + (1 - Pb1) * muPBNRPC #BBC11
expected_pref[16] <- Pb1 * muPBR + (1 - Pb1) * muPBNRPD #BBD11

### Solve for "true parameters" for PRPP-SMART DTR model -- depends on the chosen model! ###
mu_mat <- matrix(expected_pref, nrow=16)
mu_mat_trad <- matrix(expected_pref[1:4], nrow=4)
#Note: Always need to specify contrast matrix so that DTRs are in the following order
#Note: geeglm will reorder coefficients so that interactions are last, so specify formula string with main effects first
#AAC00, AAD00, BBC00, BBD00, AAC01, AAD01, BBC01, BBD01, AAC10, AAD10, BBC10, BBD10, AAC11, AAD11, BBC11, BBD11

## Model 1:
contrast_dtr_M1 <- matrix(c(1,1,0,1,0,1, #AC00
                         1,1,0,0,0,0, #AD00
                         1,0,0,1,0,0, #BC00
                         1,0,0,0,0,0, #BD00
                         1,1,0,1,1,1, #AC01
                         1,1,0,0,1,0, #AD01
                         1,0,0,1,1,0, #BC01
                         1,0,0,0,1,0, #BD01
                         1,1,1,1,0,1, #AC10
                         1,1,1,0,0,0, #AD10
                         1,0,1,1,0,0, #BC10
                         1,0,1,0,0,0, #BD10
                         1,1,1,1,1,1, #AC11
                         1,1,1,0,1,0, #AD11
                         1,0,1,1,1,0, #BC11
                         1,0,1,0,1,0 #BD11
),nrow = 16, ncol=6, byrow = TRUE)
n_param_M1 = 6
param_names_M1 = c("alpha1", "beta1", "delta1", "theta1", "delta2", "gamma1")
param_compare_M1 = c(1,2,4,6) #index of parameters to compare to traditional analysis
formula_string_M1 = "Y ~ T1 + S1_Preference + T2 + S2_Preference + T1:T2"
formula_M1 = as.formula(formula_string_M1)
true_parameters_M1 <-solve(t(contrast_dtr_M1)%*%(contrast_dtr_M1))%*%t(contrast_dtr_M1)%*%mu_mat
true_parameterDTR_mat_M1 <- matrix(c(true_parameters_M1, expected_pref), ncol = 1)
rownames(true_parameterDTR_mat_M1) <-c(param_names_M1, dtr_names)

## Model 2:
contrast_dtr_M2 <- matrix(c(1,1,1,1,0,0,0,0,0, #AC00
                         1,1,0,0,0,0,0,0,0, #AD00
                         1,0,1,0,0,0,0,0,0, #BC00
                         1,0,0,0,0,0,0,0,0, #BD00
                         1,1,1,1,0,0,1,0,0, #AC01
                         1,1,0,0,0,0,0,1,0, #AD01
                         1,0,1,0,0,0,1,0,0, #BC01
                         1,0,0,0,0,0,0,1,0, #BD01
                         1,1,1,1,1,0,0,0,0, #AC10
                         1,1,0,0,1,0,0,0,0, #AD10
                         1,0,1,0,0,1,0,0,0, #BC10
                         1,0,0,0,0,1,0,0,0, #BD10
                         1,1,1,1,1,0,1,0,1, #AC11
                         1,1,0,0,1,0,0,1,1, #AD11
                         1,0,1,0,0,1,1,0,1, #BC11
                         1,0,0,0,0,1,0,1,1 #BD11
),nrow = 16, ncol=9, byrow = TRUE)
n_param_M2 = 9
param_names_M2 = c("alpha2", "beta2", "theta2",  "gamma2", "omegaA", "omegaB", "omegaC", "omegaD", "phi")
param_compare_M2 =  c(1,2,4,6) #index of parameters to compare to traditional analysis
#Note: use T1_A, T1_B, T2_C, T2_D so that I can code parameters in the right order
#otherwise, gee will automatically put main effects first
formula_string_M2 = "Y ~ T1 + T2 + T1:T2 + S1_Preference:T1_A + S1_Preference:T1_B + S2_Preference:T2_C + S2_Preference:T2_D + S1_Preference:S2_Preference"
formula_M2 = as.formula(formula_string_M2)
true_parameters_M2 <-solve(t(contrast_dtr_M2)%*%(contrast_dtr_M2))%*%t(contrast_dtr_M2)%*%mu_mat
true_parameterDTR_mat_M2 <- matrix(c(true_parameters_M2, expected_pref), ncol = 1)
rownames(true_parameterDTR_mat_M2) <-c(param_names_M2, dtr_names)

## Traditional analysis
contrast_dtr_trad <- matrix(c(1,1,1,1, #AC00
                              1,1,0,0, #AD00
                              1,0,1,0, #BC00
                              1,0,0,0 #BD00
),nrow = 4, ncol=4, byrow = TRUE)
param_names_trad = c("alpha3", "beta3", "theta3", "gamma3")
formula_string_trad = "Y ~ T1 + T2 + T1:T2"
formula_trad = as.formula(formula_string_trad)
true_param_trad <- solve(contrast_dtr_trad)%*%mu_mat_trad
true_parameterDTR_mat_trad<- matrix(c(true_param_trad, expected_pref[1:4]), ncol = 1)
rownames(true_parameterDTR_mat_trad) <- c("alpha3", "beta3", "theta3", "gamma3", "AAC00", "AAD00", "BBC00", "BBD00")

###### SIMULATION ######
## DTR analysis results -- Model 1
DTR_hat_M1 <- matrix(NA, nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation 
param_hat_M1 <- matrix(NA, nrow=n_param_M1, ncol=n.sim) # matrix to store parameter estimates per simulation 
variance_param_hat_M1 <- matrix(NA, nrow=n_param_M1, ncol=n.sim) # matrix to store sandwhich variance of parameter estimates per simulation
variance_dtr_hat_M1 <- matrix(NA, nrow=16, ncol=n.sim) # matrix to store delta method DTR variance estimates per simulation
ci_hat_M1 <- matrix(NA, nrow=n_param_M1+16, ncol=n.sim) # matrix to store whether ci covers truth per simulation  
rownames(DTR_hat_M1) <- rownames(variance_dtr_hat_M1) <- dtr_names
rownames(param_hat_M1) <- rownames(variance_param_hat_M1) <- param_names_M1
rownames(ci_hat_M1) <- c(param_names_M1, dtr_names)
## DTR analysis results -- Model 2
DTR_hat_M2 <- matrix(NA, nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation 
param_hat_M2 <- matrix(NA, nrow=n_param_M2, ncol=n.sim) # matrix to store parameter estimates per simulation 
variance_param_hat_M2 <- matrix(NA, nrow=n_param_M2, ncol=n.sim) # matrix to store sandwhich variance of parameter estimates per simulation
variance_dtr_hat_M2 <- matrix(NA, nrow=16, ncol=n.sim) # matrix to store delta method DTR variance estimates per simulation
ci_hat_M2 <- matrix(NA, nrow=n_param_M2+16, ncol=n.sim) # matrix to store whether ci covers truth per simulation  
rownames(DTR_hat_M2) <- rownames(variance_dtr_hat_M2) <- dtr_names
rownames(param_hat_M2) <- rownames(variance_param_hat_M2) <- param_names_M2
rownames(ci_hat_M2) <- c(param_names_M2, dtr_names)
## DTR analysis results -- traditional analysis
DTR_hat_trad <- matrix(NA, nrow = 4, ncol = n.sim) # matrix to store preference DTR estimate per simulation 
param_hat_trad <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store parameter estimates per simulation 
variance_param_hat_trad <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store sandwhich variance of parameter estimates per simulation
variance_dtr_hat_trad <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store delta method DTR variance estimates per simulation
ci_hat_trad <- matrix(NA, nrow=8, ncol=n.sim) # matrix to store whether ci covers truth per simulation  
rownames(DTR_hat_trad) <- rownames(variance_dtr_hat_trad) <- dtr_names[1:4]
rownames(param_hat_trad) <- rownames(variance_param_hat_trad) <- param_names_trad
rownames(ci_hat_trad) <- c(param_names_trad, dtr_names[1:4])

# DTR & Pathway sample sizes
n.DTR <- matrix(NA, nrow = 16, ncol = n.sim) # matrix to store sample size for each DTR path per simulation
n_RAR <- rep(NA, n.sim) # pathway n: randomized A responded
n_RBR <- rep(NA, n.sim) # pathway n: randomized B responded
n_RANRRC <- rep(NA, n.sim) # pathway n: randomized A no response randomized C
n_RANRRD <- rep(NA, n.sim) # pathway n: randomized A no response randomized D
n_RBNRRC <- rep(NA, n.sim) # pathway n: randomized B no response randomized C
n_RBNRRD <- rep(NA, n.sim) # pathway n: randomized B no response randomized D
#01
n_RANRPC <- rep(NA, n.sim) # pathway n: randomized A no response preferred C
n_RANRPD <- rep(NA, n.sim) # pathway n: randomized A no response preferred D
n_RBNRPC <- rep(NA, n.sim) # pathway n: randomized B no response preferred C
n_RBNRPD <- rep(NA, n.sim) # pathway n: randomized B no response preferred D
#10
n_PAR <- rep(NA, n.sim) # pathway n: preferred A responded
n_PBR <- rep(NA, n.sim) # pathway n: preferred B responded
n_PANRRC <- rep(NA, n.sim) # pathway n: preferred A no response randomized C
n_PANRRD <- rep(NA, n.sim)  # pathway n: preferred A no response randomized D
n_PBNRRC <- rep(NA, n.sim) # pathway n: preferred B no response randomized C
n_PBNRRD <- rep(NA, n.sim) # pathway n: preferred B no response randomized D
#11
n_PANRPC <- rep(NA, n.sim) # pathway n: preferred A no response preferred C
n_PANRPD <- rep(NA, n.sim) # pathway n: preferred A no response preferred D
n_PBNRPC <- rep(NA, n.sim) # pathway n: preferred B no response preferred C
n_PBNRPD <- rep(NA, n.sim) # pathway n: preferred B no response preferred D

##### Run simulations #####
seed = 0 #current simulation attempt
num_skip = 0 #count number of skipped simulations
iterations = 0 #current number of simulations for analysis
i = 0
# Scenario a you need 507 attempts to get 500 simulations
# Scenario b you need 549 attempts to get 500 simulations
# Scenario c you need 718 attempts to get 500 simulations

while(iterations < n.sim){
  seed = seed + 1
  set.seed(seed+10000)
  
  #### Generate data and check if simulation needs to be skipped ####
  data <- generate_data(N=N, pNP_target=pNP_target, pTheta_target=pTheta_target, pNP2_target=pNP2_target, pTheta2_target=pTheta2_target, Pa=Pa, Pb=Pb, Pa1=Pa1, Pb1=Pb1, muPAR=muPAR, muRAR=muRAR, muPBR=muPBR, muRBR=muRBR, muPANRPC=muPANRPC, muPANRRC=muPANRRC, muPANRPD=muPANRPD, muPANRRD=muPANRRD, muRANRPC=muRANRPC, muRANRRC=muRANRRC, muRANRPD=muRANRPD, muRANRRD=muRANRRD, muPBNRPC=muPBNRPC, muPBNRRC=muPBNRRC, muPBNRPD=muPBNRPD, muPBNRRD=muPBNRRD, muRBNRPC=muRBNRPC, muRBNRRC=muRBNRRC, muRBNRPD=muRBNRPD, muRBNRRD=muRBNRRD, sigma2=sigma2)
  trialpath_df <- data[[1]] %>% dplyr::group_by(Trial_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
  if (nrow(trialpath_df) < 20) {num_skip = num_skip + 1 ; next} # check to make sure data in each pathway 
  # If the simulation is not skipped, increase number of analysis iterations
  iterations = iterations + 1
  i = iterations #index for storing data
  
  #### Create analysis ready data sets, store basic sample size info ####
  df <- data[[2]] # data for prpp_smart full analysis (all subjects) replicated and has weights
  df <-df[order(df$id),] # sort data by id for gee statement
  df$T1_A = ifelse(df$T1 == 1,1,0)
  df$T1_B = ifelse(df$T1 == 0,1,0)
  df$T2_C = ifelse(df$T2 == 1,1,0)
  df$T2_D = ifelse(df$T2 == 0,1,0)
  
  ## create data set for traditional analysis -- treatment indifferent participants only ## 
  df2 <- data[[1]] # relabeled raw data (not replicated has calculated weights)
  trad.ind <- which(df2$Trial_Path == "RAR" | df2$Trial_Path == "RANRRC" | df2$Trial_Path == "RANRRD"| df2$Trial_Path == "RBR" | df2$Trial_Path == "RBNRRC" | df2$Trial_Path == "RBNRRD") # only want tradional randomized SMART sujbects
  df_rand <- df2[trad.ind,]
  
  # replicate data - only twice for traditional analysis 
  # 1st dataset of responders setting T2=1, 
  datareps1 <- df_rand[df_rand$R==1,] 
  datareps1$T2 <- 1
  
  # 2nd dataset of responders setting T2=0, 
  datareps2 <- df_rand[df_rand$R==1,] 
  datareps2$T2 <- 0
  
  # dataset for non-responders
  datanoresp <- df_rand[df_rand$R==0,]
  
  # replicated data
  replicated_dat_trad <- rbind(datareps1,datareps2,datanoresp)
  
  # create data used in traditional analysis 
  analysis_data_trad <- replicated_dat_trad %>% 
    dplyr::select(id, X1, X2, Y, w, T1, T2, R, Trial_Path)
  
  analysis_data_trad <-analysis_data_trad[order(analysis_data_trad$id),] # sort data by id for gee statement
  
  # DTR N Track
  n.DTR[,i] = data[[3]]
  #00
  n_RAR[i] <- sum(data[[1]]$Trial_Path == "RAR")
  n_RBR[i] <- sum(data[[1]]$Trial_Path == "RBR")
  n_RANRRC[i] <- sum(data[[1]]$Trial_Path == "RANRRC")
  n_RANRRD[i] <- sum(data[[1]]$Trial_Path == "RANRRD")
  n_RBNRRC[i] <- sum(data[[1]]$Trial_Path == "RBNRRC")
  n_RBNRRD[i] <- sum(data[[1]]$Trial_Path == "RBNRRD")
  #01
  n_RANRPC[i] <- sum(data[[1]]$Trial_Path == "RANRPC")
  n_RANRPD[i] <- sum(data[[1]]$Trial_Path == "RANRPD")
  n_RBNRPC[i] <- sum(data[[1]]$Trial_Path == "RBNRPC")
  n_RBNRPD[i] <- sum(data[[1]]$Trial_Path == "RBNRPD")
  #10
  n_PAR[i] <- sum(data[[1]]$Trial_Path == "PAR")
  n_PBR[i] <- sum(data[[1]]$Trial_Path == "PBR")
  n_PANRRC[i] <- sum(data[[1]]$Trial_Path == "PANRRC")
  n_PANRRD[i] <- sum(data[[1]]$Trial_Path == "PANRRD")
  n_PBNRRC[i] <- sum(data[[1]]$Trial_Path == "PBNRRC")
  n_PBNRRD[i] <- sum(data[[1]]$Trial_Path == "PBNRRD")
  #11
  n_PANRPC[i] <- sum(data[[1]]$Trial_Path == "PANRPC")
  n_PANRPD[i] <- sum(data[[1]]$Trial_Path == "PANRPD")
  n_PBNRPC[i] <- sum(data[[1]]$Trial_Path == "PBNRPC")
  n_PBNRPD[i] <- sum(data[[1]]$Trial_Path == "PBNRPD")
  
  #### DTR Analyses ####
  
  ### PRPP-SMART DTR Analysis, Model 1 ###
  # use model formula loaded from model file
  gee.fit.1 <- geeglm(formula_M1, id=id, weights=w, family = gaussian, corstr = "independence", data = df)
  
  # Parameter & DTR estimates
  param_hat_M1[,i] <- gee.fit.1$coefficients
  DTR_hat_M1[,i] = contrast_dtr_M1 %*% param_hat_M1[,i]
  
  # Robust variances of each parameter 
  variance_param_hat_M1[,i] <- diag(gee.fit.1$geese$vbeta)
  
  # DTR variances (sum of parameters, so no longer need delta method)
  variance_dtr_hat_M1[,i] = apply(contrast_dtr_M1, 1, function(x) t(as.matrix(x)) %*% vcov(gee.fit.1) %*% as.matrix(x))
  
  # Parameter & DTR CI coverage
  ci_mat <- confint.geeglm(gee.fit.1) # Parameter coverage
  lwr_dtr <- c()
  upper_dtr <- c()
  for(d in 1:16){
    lwr_dtr[d] <- DTR_hat_M1[d,i] - qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_M1[d,i])
    upper_dtr[d] <- DTR_hat_M1[d,i] + qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_M1[d,i])
  }
  dtr_ci_mat <- cbind(lwr_dtr, upper_dtr)
  ci_mat <- rbind(ci_mat, dtr_ci_mat)
  param_in_ci <- data.table::between(true_parameterDTR_mat_M1[,1], ci_mat[, 1], ci_mat[, 2])
  ci_hat_M1[,i] <- param_in_ci
  
  ### PRPP-SMART DTR Analysis, Model 2 ###
  # use model formula loaded from model file
  gee.fit.2 <- geeglm(formula_M2, id=id, weights=w, family = gaussian, corstr = "independence", data = df)
  
  # Parameter & DTR estimates
  param_hat_M2[,i] <- gee.fit.2$coefficients
  DTR_hat_M2[,i] = contrast_dtr_M2 %*% param_hat_M2[,i]
  
  # Robust variances of each parameter 
  variance_param_hat_M2[,i] <- diag(gee.fit.2$geese$vbeta)
  
  # DTR variances (sum of parameters, so no longer need delta method)
  variance_dtr_hat_M2[,i] = apply(contrast_dtr_M2, 1, function(x) t(as.matrix(x)) %*% vcov(gee.fit.2) %*% as.matrix(x))
  
  # Parameter & DTR CI coverage
  ci_mat <- confint.geeglm(gee.fit.2) # Parameter coverage
  lwr_dtr <- c()
  upper_dtr <- c()
  for(d in 1:16){
    lwr_dtr[d] <- DTR_hat_M2[d,i] - qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_M2[d,i])
    upper_dtr[d] <- DTR_hat_M2[d,i] + qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_M2[d,i])
  }
  dtr_ci_mat <- cbind(lwr_dtr, upper_dtr)
  ci_mat <- rbind(ci_mat, dtr_ci_mat)
  param_in_ci <- data.table::between(true_parameterDTR_mat_M2[,1], ci_mat[, 1], ci_mat[, 2])
  ci_hat_M2[,i] <- param_in_ci
  
  ### Traditional PRPP-SMART ANALYSIS (randomized randomized subjects only) ###
  
  # fit model
  gee.fit.t <- geeglm(formula_trad, id=id, weights=w, family = gaussian, corstr = "independence", data = analysis_data_trad)
  
  # Parameter & DTR estimates
  param_hat_trad[,i] <- gee.fit.t$coefficients
  DTR_hat_trad[,i] = contrast_dtr_trad %*% param_hat_trad[,i]
  
  # Robust variances of each parameter 
  variance_param_hat_trad[,i] <- diag(gee.fit.t$geese$vbeta)
  
  ## DTR variances (no delta method needed)
  variance_dtr_hat_trad[,i] = apply(contrast_dtr_trad, 1, function(x) t(as.matrix(x)) %*% vcov(gee.fit.t) %*% as.matrix(x))
  
  # CI coverage 
  ci_mat <- confint.geeglm(gee.fit.t) # parameter CIs
  lwr_dtr <- c()
  upper_dtr <- c()
  for(d in 1:4){
    lwr_dtr[d] <- DTR_hat_trad[d,i] - qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_trad[d,i])
    upper_dtr[d] <- DTR_hat_trad[d,i] + qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_trad[d,i])
  }
  dtr_ci_mat <- cbind(lwr_dtr, upper_dtr)
  ci_mat <- rbind(ci_mat, dtr_ci_mat)
  param_in_ci <- data.table::between(true_parameterDTR_mat_trad[,1], ci_mat[, 1], ci_mat[, 2])
  ci_hat_trad[,i] <- param_in_ci
  
  # PRPP-SMART Main Effects Analysis -- see https://github.com/mariwank/PRPP-SMART
}

# number simulations skipped
num_skip_total <- num_skip

###### Evaluation & File Output ######
### Model 1
## Parameters
param_hat_avg_M1 = c()
param_sd_hat_M1 = c()
param_avg_sd_M1 = c()

for(i in 1:length(param_names_M1)){
  param_hat_avg_M1[i] <- round(mean(param_hat_M1[i,], na.rm = TRUE),4)
  param_sd_hat_M1[i] <- round(sd(param_hat_M1[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 parameter estimates)
  param_avg_sd_M1[i] <- round(sqrt(mean(variance_param_hat_M1[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(parameter hat))) equivalent to mean(sqrt(var(parameter hat)))
}
# calculate bias
bias_param_hat_avg_M1 <- round(param_hat_avg_M1-true_parameters_M1,4)
# calculate rMSE
rMSE_param_M1 <- sqrt(param_sd_hat_M1^2 + bias_param_hat_avg_M1^2)

param_results_tbl_M1 <- data.frame(Parameter=param_names_M1,True_Parameter=round(true_parameters_M1,4), 
                                Param_Hat_Avg=param_hat_avg_M1, 
                                Bias = bias_param_hat_avg_M1, 
                                SE = param_sd_hat_M1, # sd of 500 parameter estimates
                                Avg_se = param_avg_sd_M1, # mean of 500 parameter variances then square root
                                rMSE = round(rMSE_param_M1,4)) # rMSE using se calculated from 500 parameter estimates 
write.csv(param_results_tbl_M1, paste0(outdir, "/ParamResults_Model1_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)


## DTRs
# calculate bias
DTR_hat_avg_M1 = c()
DTR_sd_hat_M1 = c()
DTR_avg_sd_M1 = c()
DTR_avg_n = c()

for(i in 1:16){
  DTR_hat_avg_M1[i] <- round(mean(DTR_hat_M1[i,], na.rm = TRUE),4)
  DTR_sd_hat_M1[i] <- round(sd(DTR_hat_M1[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 DTR estimates)
  DTR_avg_sd_M1[i] <- round(sqrt(mean(variance_dtr_hat_M1[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(dtr hat))) equivalent to mean(sqrt(var(dtr hat)))
  DTR_avg_n[i] <- round(mean(n.DTR[i,], na.rm = TRUE),1)
}
# calculate bias
DTR_bias_M1 = round(DTR_hat_avg_M1 - expected_pref,4)
# calculate rMSE
rMSE_DTR_M1 <- sqrt(DTR_sd_hat_M1^2 + DTR_bias_M1^2)

DTR_results_tbl_M1 <- data.frame(DTR=dtr_names,True_DTR=expected_pref,
                              DTR_Hat_Avg=DTR_hat_avg_M1, 
                              Bias=DTR_bias_M1, 
                              SE = DTR_sd_hat_M1, # sd of 500 DTR estimates
                              Avg_se = DTR_avg_sd_M1, # mean of 500 DTR variances then square root
                              rMSE = round(rMSE_DTR_M1,4), # rMSE using sd calculated from 500 DTR estimates   
                              DTR_AvgN = DTR_avg_n) 
write.csv(DTR_results_tbl_M1, paste0(outdir, "/DTRResults_Model1_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

## Coverage Rates
ci_coverage_M1 = c()
for(i in 1:nrow(ci_hat_M1)){
  ci_coverage_M1[i] <- mean(ci_hat_M1[i,], na.rm = TRUE)
}
ci_results_tbl_M1 <- data.frame(Parameter_DTR=c(param_names_M1, dtr_names), 
                             CI_Coverage=round(ci_coverage_M1,3))
write.csv(ci_results_tbl_M1, paste0(outdir, "/CIResults_Model1_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

### Model 2
## Parameters
param_hat_avg_M2 = c()
param_sd_hat_M2 = c()
param_avg_sd_M2 = c()

for(i in 1:length(param_names_M2)){
  param_hat_avg_M2[i] <- round(mean(param_hat_M2[i,], na.rm = TRUE),4)
  param_sd_hat_M2[i] <- round(sd(param_hat_M2[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 parameter estimates)
  param_avg_sd_M2[i] <- round(sqrt(mean(variance_param_hat_M2[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(parameter hat))) equivalent to mean(sqrt(var(parameter hat)))
}
# calculate bias
bias_param_hat_avg_M2 <- round(param_hat_avg_M2-true_parameters_M2,4)
# calculate rMSE
rMSE_param_M2 <- sqrt(param_sd_hat_M2^2 + bias_param_hat_avg_M2^2)

param_results_tbl_M2 <- data.frame(Parameter=param_names_M2,True_Parameter=round(true_parameters_M2,4), 
                                   Param_Hat_Avg=param_hat_avg_M2, 
                                   Bias = bias_param_hat_avg_M2, 
                                   SE = param_sd_hat_M2, # sd of 500 parameter estimates
                                   Avg_se = param_avg_sd_M2, # mean of 500 parameter variances then square root
                                   rMSE = round(rMSE_param_M2,4)) # rMSE using se calculated from 500 parameter estimates 
write.csv(param_results_tbl_M2, paste0(outdir, "/ParamResults_Model2_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)


## DTRs
# calculate bias
DTR_hat_avg_M2 = c()
DTR_sd_hat_M2 = c()
DTR_avg_sd_M2 = c()
DTR_avg_n = c()

for(i in 1:16){
  DTR_hat_avg_M2[i] <- round(mean(DTR_hat_M2[i,], na.rm = TRUE),4)
  DTR_sd_hat_M2[i] <- round(sd(DTR_hat_M2[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 DTR estimates)
  DTR_avg_sd_M2[i] <- round(sqrt(mean(variance_dtr_hat_M2[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(dtr hat))) equivalent to mean(sqrt(var(dtr hat)))
  DTR_avg_n[i] <- round(mean(n.DTR[i,], na.rm = TRUE),1)
}
# calculate bias
DTR_bias_M2 = round(DTR_hat_avg_M2 - expected_pref,4)
# calculate rMSE
rMSE_DTR_M2 <- sqrt(DTR_sd_hat_M2^2 + DTR_bias_M2^2)

DTR_results_tbl_M2 <- data.frame(DTR=dtr_names,True_DTR=expected_pref,
                                 DTR_Hat_Avg=DTR_hat_avg_M2, 
                                 Bias=DTR_bias_M2, 
                                 SE = DTR_sd_hat_M2, # sd of 500 DTR estimates
                                 Avg_se = DTR_avg_sd_M2, # mean of 500 DTR variances then square root
                                 rMSE = round(rMSE_DTR_M2,4), # rMSE using sd calculated from 500 DTR estimates   
                                 DTR_AvgN = DTR_avg_n) 
write.csv(DTR_results_tbl_M2, paste0(outdir, "/DTRResults_Model2_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

## Coverage Rates
ci_coverage_M2 = c()
for(i in 1:nrow(ci_hat_M2)){
  ci_coverage_M2[i] <- mean(ci_hat_M2[i,], na.rm = TRUE)
}
ci_results_tbl_M2 <- data.frame(Parameter_DTR=c(param_names_M2, dtr_names), 
                                CI_Coverage=round(ci_coverage_M2,3))
write.csv(ci_results_tbl_M2, paste0(outdir, "/CIResults_Model2_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

#### Traditional DTR Analysis ####
## Parameters
param_hat_avg_trad = c()
param_sd_hat_trad = c()
param_avg_sd_trad = c()

for(i in 1:length(param_names_trad)){
  param_hat_avg_trad[i] <- round(mean(param_hat_trad[i,], na.rm = TRUE),4)
  param_sd_hat_trad[i] <- round(sd(param_hat_trad[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 parameter estimates)
  param_avg_sd_trad[i] <- round(sqrt(mean(variance_param_hat_trad[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(parameter hat))) equivalent to mean(sqrt(var(parameter hat)))
}

# calculate bias
bias_param_hat_avg_trad <- round(param_hat_avg_trad-true_param_trad,4)
# calculate rMSE 
rMSE_param_trad <- sqrt(param_sd_hat_trad^2 + bias_param_hat_avg_trad^2)

param_results_tbl_trad <- data.frame(Parameter=param_names_trad,True_Parameter=round(true_param_trad,4), 
                                  Param_Hat_Avg=param_hat_avg_trad, 
                                  Bias = bias_param_hat_avg_trad, 
                                  SE = param_sd_hat_trad, # sd of 500 parameter estimates
                                  Avg_se = param_avg_sd_trad, # mean of 500 parameter variances then square root
                                  rMSE = round(rMSE_param_trad,4))# rMSE using se calculated from 500 parameter estimates
write.csv(param_results_tbl_trad, paste0(outdir, "/ParamResults_Trad_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

## DTR 
DTR_hat_avg_trad = c()
DTR_sd_hat_trad = c()
DTR_avg_sd_trad = c()

for(i in 1:4){
  DTR_hat_avg_trad[i] <- round(mean(DTR_hat_trad[i,], na.rm = TRUE),4)
  DTR_sd_hat_trad[i] <- round(sd(DTR_hat_trad[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 DTR estimates)
  DTR_avg_sd_trad[i] <- round(sqrt(mean(variance_dtr_hat_trad[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(dtr hat))) equivalent to mean(sqrt(var(dtr hat)))
}

# calculate bias
DTR_bias_trad = round(DTR_hat_avg_trad - expected_pref[1:4],4)
# calculate rMSE
rMSE_DTR_trad <- sqrt(DTR_sd_hat_trad^2 + DTR_bias_trad^2)

DTR_results_tbl_trad <- data.frame(DTR=dtr_names[1:4],True_DTR=expected_pref[1:4],
                                DTR_Hat_Avg=DTR_hat_avg_trad, 
                                Bias=DTR_bias_trad, 
                                SE = DTR_sd_hat_trad, # sd of 500 DTR estimates
                                Avg_se = DTR_avg_sd_trad, # mean of 500 DTR variances then square root
                                rMSE = round(rMSE_DTR_trad,4), # rMSE using sd calculated from 500 DTR estimates
                                DTR_AvgN = DTR_avg_n[1:4]) 
write.csv(DTR_results_tbl_trad, paste0(outdir, "/DTRResults_Trad_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

## Coverage Rates
ci_coverage_trad = c()

for(i in 1:nrow(ci_hat_trad)){
  ci_coverage_trad[i] <- mean(ci_hat_trad[i,], na.rm = TRUE)
}
ci_results_tbl_trad <- data.frame(Parameter_DTR=c(param_names_trad,dtr_names[1:4]), 
                               CI_Coverage=round(ci_coverage_trad,3))
write.csv(ci_results_tbl_trad, paste0(outdir, "/CIResults_Trad_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)

## Pathway N
#00
n_RAR_avg <- mean(n_RAR, na.rm=T)
n_RBR_avg <- mean(n_RBR, na.rm=T)
n_RANRRC_avg <- mean(n_RANRRC, na.rm=T)
n_RANRRD_avg <- mean(n_RANRRD, na.rm=T)
n_RBNRRC_avg <- mean(n_RBNRRC, na.rm=T)
n_RBNRRD_avg <- mean(n_RBNRRD, na.rm=T)
#01
n_RANRPC_avg <- mean(n_RANRPC, na.rm=T)
n_RANRPD_avg <- mean(n_RANRPD, na.rm=T)
n_RBNRPC_avg <- mean(n_RBNRPC, na.rm=T)
n_RBNRPD_avg <- mean(n_RBNRPD, na.rm=T)
#10
n_PAR_avg <- mean(n_PAR, na.rm=T)
n_PBR_avg <- mean(n_PBR, na.rm=T)
n_PANRRC_avg <- mean(n_PANRRC, na.rm=T)
n_PANRRD_avg <- mean(n_PANRRD, na.rm=T)
n_PBNRRC_avg <- mean(n_PBNRRC, na.rm=T)
n_PBNRRD_avg <- mean(n_PBNRRD, na.rm=T)
#11
n_PANRPC_avg <- mean(n_PANRPC, na.rm=T)
n_PANRPD_avg <- mean(n_PANRPD, na.rm=T)
n_PBNRPC_avg <- mean(n_PBNRPC, na.rm=T)
n_PBNRPD_avg <- mean(n_PBNRPD, na.rm=T)

pathway_n_dtr_df <- data.frame(Pathway = c("RAR",
                                           "RBR",
                                           "PAR",
                                           "PBR",
                                           "RANRRC", 
                                           "RANRRD", 
                                           "RBNRRC",
                                           "RBNRRD",
                                           "RANRPC", 
                                           "RANRPD", 
                                           "RBNRPC",
                                           "RBNRPD",
                                           "PANRRC", 
                                           "PANRRD", 
                                           "PBNRRC",
                                           "PBNRRD",
                                           "PANRPC", 
                                           "PANRPD", 
                                           "PBNRPC",
                                           "PBNRPD"),
                               Avg_N = c(n_RAR_avg,
                                         n_RBR_avg,
                                         n_PAR_avg,
                                         n_PBR_avg,
                                         n_RANRRC_avg,
                                         n_RANRRD_avg,
                                         n_RBNRRC_avg,
                                         n_RBNRRD_avg,
                                         n_RANRPC_avg,
                                         n_RANRPD_avg,
                                         n_RBNRPC_avg,
                                         n_RBNRPD_avg,
                                         n_PANRRC_avg,
                                         n_PANRRD_avg,
                                         n_PBNRRC_avg,
                                         n_PBNRRD_avg,
                                         n_PANRPC_avg,
                                         n_PANRPD_avg,
                                         n_PBNRPC_avg,
                                         n_PBNRPD_avg))

write.csv(pathway_n_dtr_df, paste0(outdir, "/PathwayN_Scenario", type, scenario, "_", size, "Effects_Sigma2_", sigma2, "_N_", N, "_NSim", n.sim, ".csv"), row.names=FALSE)
