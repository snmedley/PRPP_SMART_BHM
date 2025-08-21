# Install JAGS before rjags!
library(rjags)
library(coda) # For MCMC convergence diagnostics

# Choose home directory where R code is stored. 
# Data generation files must be stored in a subdirectory named "DataGeneration"
# Files with values of prior parameters and BUGS files must be stored in a subdirectory named "JAGSFiles"
homedir = "C:/Users/snmed/OneDrive/Documents/GSRA/Project_PCORI_PRPP/DynamicBorrowing/GitHubCode"
# Choose a directory for output files. The directory is set to homedir by default. 
outdir = homedir

##### Choose Data Generation Settings #####
# Choose sample size
N = 500
# Load data generation function (depends on chosen model)
setwd(paste0(homedir, "/DataGeneration"))
source("PRPP_SMART_DataGen_cts.R")
# Choose sigma2 (data variability at trial pathway level) 
sigma2 = 36
# Choose preference rate scenario
scenario = "a" # "a", "b", "c", "d", "e", or "f"
load(paste0("scenario_", scenario, ".RData"))
# Choose preference augmented DTR effect type and effect size (delta)
type = 1 # 1, 2, 3, or 4
size = "small" # "small", "moderate", or "large"
# Note: all combinations are valid except type 4 with a small effect size
load(paste0("type", type, "_", size, ".RData"))

##### Calculate True Parameter and DTR Values #####
#### DTRs ####
# True DTR values for PRPP-SMART -- depends on Pa, Pa1, Pb, Pb1 from scenario file
# Also depends on pathway means specified in preference setting file
# Does NOT depend on the chosen model
dtr.names = c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
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

##### Choose Prior Settings and MCMC Settings #####
setwd("../JAGSFiles")
load("Mu_Prior1.RData")
tau_sd = 3/2 # 3/2 for a prior assuming substantial heterogeneity, 3 for large heterogeneity
## JAGS settings
# Note: JAGS will output parameters in alphabetical order in model summary
pars_trad = c("theta_A_0", "theta_A_1", "theta_B_0", "theta_B_1",
              "theta_AC_00", "theta_AC_01", "theta_AC_10", "theta_AC_11",   
              "theta_AD_00", "theta_AD_01", "theta_AD_10", "theta_AD_11",  
              "theta_BC_00", "theta_BC_01", "theta_BC_10", "theta_BC_11",  
              "theta_BD_00", "theta_BD_01", "theta_BD_10", "theta_BD_11",  
              "theta_AAC00", "theta_AAD00", "theta_BBC00", "theta_BBD00", # parameters to monitor
              "theta_AAC01", "theta_AAD01", "theta_BBC01",  "theta_BBD01",
              "theta_AAC10", "theta_AAD10", "theta_BBC10", "theta_BBD10",
              "theta_AAC11", "theta_AAD11", "theta_BBC11", "theta_BBD11")
pars_full = c("p_A_0", "p_A_1", "p_B_0", "p_B_1", "tau_NR", "tau_R", 
              "theta_A_0", "theta_A_1", "theta_B_0", "theta_B_1",
              "theta_AC_00", "theta_AC_01", "theta_AC_10", "theta_AC_11",   
              "theta_AD_00", "theta_AD_01", "theta_AD_10", "theta_AD_11",  
              "theta_BC_00", "theta_BC_01", "theta_BC_10", "theta_BC_11",  
              "theta_BD_00", "theta_BD_01", "theta_BD_10", "theta_BD_11",  
              "theta_AAC00", "theta_AAD00", "theta_BBC00", "theta_BBD00", # parameters to monitor
              "theta_AAC01", "theta_AAD01", "theta_BBC01",  "theta_BBD01",
              "theta_AAC10", "theta_AAD10", "theta_BBC10", "theta_BBD10",
              "theta_AAC11", "theta_AAD11", "theta_BBC11", "theta_BBD11")
n.adapt <- 5000 
n.burnin <- 5000
MCMC_SAMPLE <- 100000 #number of posterior samples AFTER adaptation and burn-in
n.thin <- 50 #chains should be 2000 draws long each
n_MCMC_chain <- 3

###### SIMULATION ######
seed = 1234 
set.seed(seed)
data <- generate_data(N=N, pNP_target=pNP_target, pTheta_target=pTheta_target, pNP2_target=pNP2_target, pTheta2_target=pTheta2_target, Pa=Pa, Pb=Pb, Pa1=Pa1, Pb1=Pb1, muPAR=muPAR, muRAR=muRAR, muPBR=muPBR, muRBR=muRBR, muPANRPC=muPANRPC, muPANRRC=muPANRRC, muPANRPD=muPANRPD, muPANRRD=muPANRRD, muRANRPC=muRANRPC, muRANRRC=muRANRRC, muRANRPD=muRANRPD, muRANRRD=muRANRRD, muPBNRPC=muPBNRPC, muPBNRRC=muPBNRRC, muPBNRPD=muPBNRPD, muPBNRRD=muPBNRRD, muRBNRPC=muRBNRPC, muRBNRRC=muRBNRRC, muRBNRPD=muRBNRPD, muRBNRRD=muRBNRRD, sigma2=sigma2)
trialpath_df <- data[[1]] %>% dplyr::group_by(Trial_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
skip = ifelse(nrow(trialpath_df) < 20, 1, 0) # skip simulation if there are positivity issues
n.DTR = data[[3]] # number of participants consistent with each DTR
df <- data[[1]] # data for prpp_smart full analysis (all subjects)
## Estimate variability in the outcome in each trial pathway
df$s = NA
df[df$Trial_Path == "RAR",]$s = sd(df[df$Trial_Path == "RAR",]$Y)
df[df$Trial_Path == "PAR",]$s = sd(df[df$Trial_Path == "PAR",]$Y)
df[df$Trial_Path == "RANRRC",]$s = sd(df[df$Trial_Path == "RANRRC",]$Y)
df[df$Trial_Path == "RANRPC",]$s = sd(df[df$Trial_Path == "RANRPC",]$Y)
df[df$Trial_Path == "PANRRC",]$s = sd(df[df$Trial_Path == "PANRRC",]$Y)
df[df$Trial_Path == "PANRPC",]$s = sd(df[df$Trial_Path == "PANRPC",]$Y)
df[df$Trial_Path == "RANRRD",]$s = sd(df[df$Trial_Path == "RANRRD",]$Y)
df[df$Trial_Path == "RANRPD",]$s = sd(df[df$Trial_Path == "RANRPD",]$Y)
df[df$Trial_Path == "PANRRD",]$s = sd(df[df$Trial_Path == "PANRRD",]$Y)
df[df$Trial_Path == "PANRPD",]$s = sd(df[df$Trial_Path == "PANRPD",]$Y)
df[df$Trial_Path == "RBR",]$s = sd(df[df$Trial_Path == "RBR",]$Y)
df[df$Trial_Path == "PBR",]$s = sd(df[df$Trial_Path == "PBR",]$Y)
df[df$Trial_Path == "RBNRRC",]$s = sd(df[df$Trial_Path == "RBNRRC",]$Y)
df[df$Trial_Path == "RBNRPC",]$s = sd(df[df$Trial_Path == "RBNRPC",]$Y)
df[df$Trial_Path == "PBNRRC",]$s = sd(df[df$Trial_Path == "PBNRRC",]$Y)
df[df$Trial_Path == "PBNRPC",]$s = sd(df[df$Trial_Path == "PBNRPC",]$Y)
df[df$Trial_Path == "RBNRRD",]$s = sd(df[df$Trial_Path == "RBNRRD",]$Y)
df[df$Trial_Path == "RBNRPD",]$s = sd(df[df$Trial_Path == "RBNRPD",]$Y)
df[df$Trial_Path == "PBNRRD",]$s = sd(df[df$Trial_Path == "PBNRRD",]$Y)
df[df$Trial_Path == "PBNRPD",]$s = sd(df[df$Trial_Path == "PBNRPD",]$Y)
## Set up data for traditional analysis (tBM) -- no pooling estimates of indifference DTRs
r_rand_A <- df[df$T1 == 1 & df$S1_Preference == 0,]$R
r_rand_B <- df[df$T1 == 0 & df$S1_Preference == 0,]$R
r_pref_A <- df[df$T1 == 1 & df$S1_Preference == 1,]$R
r_pref_B <- df[df$T1 == 0 & df$S1_Preference == 1,]$R
data_trad = c(list(length(r_rand_A), r_rand_A, length(r_rand_B), r_rand_B, 
                   length(r_pref_A), r_pref_A, length(r_pref_B), r_pref_B, 
                   nrow(df), df$T1, df$S1_Preference, df$R, df$T2, df$S2_Preference, df$Y, df$s,
                   epsilon_p, mu_A_p, mu_B_p, mu_AC_p, mu_AD_p, mu_BC_p, mu_BD_p))
names(data_trad) = c("N_rand_A", "r_rand_A", "N_rand_B", "r_rand_B",
                     "N_pref_A", "r_pref_A", "N_pref_B", "r_pref_B",
                     "N", "T1", "P1", "R", "T2", "P2", "Y", "s",
                     "epsilon_p", "mu_A_p", "mu_B_p", "mu_AC_p", "mu_AD_p", "mu_BC_p", "mu_BD_p")
## Run tBM in JAGS & collect draws from posterior
jag_trad <- rjags::jags.model(
  file = paste0("tBM.bugs"),
  data = data_trad,
  n.chains = n_MCMC_chain, n.adapt = n.adapt
)
posterior_trad <- rjags::coda.samples(  #multiple chains as mcmc.list object
  model = jag_trad,
  variable.names = pars_trad,
  n.iter = MCMC_SAMPLE,
  thin = n.thin
)
draws_trad = as.data.frame( #single set of draws from all chains. Inference is based on all draws assuming convergence
  rbind(posterior_trad[[1]], posterior_trad[[2]], posterior_trad[[3]]))
#skip convergence diagnostics for traditional analysis
## Run BHM in JAGS & collect draws from posterior
data_full = c(list(nrow(df), df$T1, df$S1_Preference, df$R, df$T2, df$S2_Preference, df$Y, df$s,
                   mu_A_p, mu_B_p, mu_AC_p, mu_AD_p, mu_BC_p, mu_BD_p,
                   omega_p, delta1_p, delta2_p,
                   epsilon_p, tau_sd))
names(data_full) = c("N", "T1", "P1", "R", "T2", "P2", "Y", "s", 
                     "mu_A_p", "mu_B_p", "mu_AC_p", "mu_AD_p", "mu_BC_p", "mu_BD_p", 
                     "omega_p", "delta1_p", "delta2_p",
                     "epsilon_p", "tau_sd")
jag <- rjags::jags.model(
  file = paste0("BHM.bugs"),
  data = data_full,
  n.chains = n_MCMC_chain, n.adapt = n.adapt
)
update(jag, n.burnin) #may need to add burn-in phase after adaptation phase, but not always necessary
posterior_full <- rjags::coda.samples(  #multiple chains as mcmc.list object
  model = jag,
  variable.names = pars_full,
  n.iter = MCMC_SAMPLE,
  thin = n.thin
)
draws_full = as.data.frame( #single set of draws from all chains. Inference is based on all draws assuming convergence
  rbind(posterior_full[[1]], posterior_full[[2]], posterior_full[[3]]))

###### RESULTS ######
setwd(outdir)
## tBM Output Tables
param_t = apply(draws_trad[,c(9:18,27:36)], 2, mean)
param_var_t = apply(draws_trad[,c(9:18,27:36)], 2, var)
param_names = c("mu_AC_00", "mu_AC_01", "mu_AC_10", "mu_AC_11", "mu_AD_00", "mu_AD_01", "mu_AD_10", "mu_AD_11", "mu_A_0", "mu_A_1",
                "mu_BC_00", "mu_BC_01", "mu_BC_10", "mu_BC_11", "mu_BD_00", "mu_BD_01", "mu_BD_10", "mu_BD_11", "mu_B_0", "mu_B_1")
# The parameter estimates in this case are the pathway means
param_tbl_t = data.frame(param_names, param_t, sqrt(param_var_t), row.names = NULL)
colnames(param_tbl_t) = c("Parameter", "Estimate", "Standard Error")
write.csv(param_tbl_t, paste0("tBM_param_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)
DTR_t = apply(draws_trad[,c(1:8,19:26)], 2, mean)
DTR_var_t = apply(draws_trad[,c(1:8,19:26)], 2, var)
DTR_quantiles_t = cbind(apply(draws_trad[,c(1:8,19:26)], 2, function(x) quantile(x, probs = 0.025)), 
                        apply(draws_trad[,c(1:8,19:26)], 2, function(x) quantile(x, probs = 0.975)))
DTR_names = c("AAC00", "AAC01", "AAC10", "AAC11", "AAC00", "AAC01", "AAC10", "AAC11",
              "BBC00", "BBC01", "BBC10", "BBC11", "BBC00", "BBC01", "BBC10", "BBC11")
DTR_tbl_t = data.frame(DTR_names, DTR_t, sqrt(DTR_var_t), DTR_quantiles_t)
colnames(DTR_tbl_t) = c("DTR", "Estimate", "Standard Error", "95% CI - Lower Bound", "95% CI - Upper Bound")
write.csv(DTR_tbl_t, paste0("tBM_DTR_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)
## BHM Output Tables
#Rhat -- have the chains converged? Rhat < 1.1 suggests convergence
R_hat = coda::gelman.diag(posterior_full)$psrf[,1] #Rhat
#neff -- measure of autocorrelation in chains. One recommendation is neff > 100
n_eff = coda::effectiveSize(posterior_full) #n.eff, summed across chains
Quantity = c("p_A_0", "p_A_1", "p_B_0", "p_B_1", "tau_NR", "tau_R", "AAC00", "AAC01", "AAC10", "AAC11", "AAC00", "AAC01", "AAC10", "AAC11", "theta_AC_00", "theta_AC_01", "theta_AC_10", "theta_AC_11", "theta_AD_00", "theta_AD_01", "theta_AD_10", "theta_AD_11", "theta_A_0", "theta_A_1",
             "BBC00", "BBC01", "BBC10", "BBC11", "BBC00", "BBC01", "BBC10", "BBC11", "theta_BC_00", "theta_BC_01", "theta_BC_10", "theta_BC_11", "theta_BD_00", "theta_BD_01", "theta_BD_10", "theta_BD_11", "theta_B_0", "theta_B_1")
diag_tbl = data.frame(Quantity, R_hat, n_eff, row.names = NULL)
write.csv(diag_tbl, paste0("BHM_diag_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)
param_hat = apply(draws_full[,c(1:6,15:24,33:42)], 2, mean)
param_var = apply(draws_full[,c(1:6,15:24,33:42)], 2, var)
param_names = c("p_A_0", "p_A_1", "p_B_0", "p_B_1", "tau_NR", "tau_R", "theta_AC_00", "theta_AC_01", "theta_AC_10", "theta_AC_11", "theta_AD_00", "theta_AD_01", "theta_AD_10", "theta_AD_11", "theta_A_0", "theta_A_1",
                "theta_BC_00", "theta_BC_01", "theta_BC_10", "theta_BC_11", "theta_BD_00", "theta_BD_01", "theta_BD_10", "theta_BD_11", "theta_B_0", "theta_B_1")
# The parameter estimates in this case are the pathway level random effects -- not the estimates of the pathway means although those can easily be derived
param_tbl = data.frame(param_names, param_hat, sqrt(param_var), row.names = NULL)
colnames(param_tbl) = c("Parameter", "Estimate", "Standard Error")
write.csv(param_tbl, paste0("BHM_param_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)
DTR_hat = apply(draws_full[,c(7:14,25:32)], 2, mean)
DTR_var = apply(draws_full[,c(7:14,25:32)], 2, var)
DTR_quantiles = cbind(apply(draws_full[,c(7:14,25:32)], 2, function(x) quantile(x, probs = 0.025)), 
                      apply(draws_full[,c(7:14,25:32)], 2, function(x) quantile(x, probs = 0.975)))
DTR_tbl = data.frame(DTR_names, DTR_hat, sqrt(DTR_var), DTR_quantiles)
colnames(DTR_tbl) = c("DTR", "Estimate", "Standard Error", "95% CI - Lower Bound", "95% CI - Upper Bound")
write.csv(DTR_tbl, paste0("BHM_DTR_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)