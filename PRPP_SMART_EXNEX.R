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
load("rMAC_Prior1.RData")
tau_sd = 3/2 # 3/2 for a prior assuming substantial heterogeneity, 3 for large heterogeneity
w = 0.5 # prior likelihood of exchangeability vs. non-exchangeability, set to be equal for all pathways
## JAGS settings
# Note: JAGS will output parameters in alphabetical order in model summary
pars_jags = c("Z", "delta1", "delta2", "omega", "tau_NR", "tau_R", 
              "theta_AAC00", "theta_AAD00", "theta_BBC00", "theta_BBD00", # parameters to monitor
              "theta_AAC01", "theta_AAD01", "theta_BBC01",  "theta_BBD01",
              "theta_AAC10", "theta_AAD10", "theta_BBC10", "theta_BBD10",
              "theta_AAC11", "theta_AAD11", "theta_BBC11", "theta_BBD11")
#use Z to compute updated mixture weights
n.adapt <- 10000 
n.burnin <- 10000
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
data_jags = c(list(nrow(df), df$T1, df$S1_Preference, df$R, df$T2, df$S2_Preference, df$Y, df$s,
                   mu_A_p, mu_B_p, mu_AC_p, mu_AD_p, mu_BC_p, mu_BD_p,
                   omega_p, delta1_p, delta2_p, epsilon_p, tau_sd, w, nu,
                   m_A_0, m_A_1, m_AC_00, m_AC_01, m_AC_10, m_AC_11, m_AD_00, m_AD_01, m_AD_10, m_AD_11,
                   m_B_0, m_B_1, m_BC_00, m_BC_01, m_BC_10, m_BC_11, m_BD_00, m_BD_01, m_BD_10, m_BD_11))
names(data_jags) = c("N", "T1", "P1", "R", "T2", "P2", "Y", "s", 
                     "mu_A_p", "mu_B_p", "mu_AC_p", "mu_AD_p", "mu_BC_p", "mu_BD_p", 
                     "omega_p", "delta1_p", "delta2_p", "epsilon_p", "tau_sd", "w", "nu",
                     "m_A_0", "m_A_1", "m_AC_00", "m_AC_01", "m_AC_10", "m_AC_11", "m_AD_00", "m_AD_01", "m_AD_10", "m_AD_11",
                     "m_B_0", "m_B_1", "m_BC_00", "m_BC_01", "m_BC_10", "m_BC_11", "m_BD_00", "m_BD_01", "m_BD_10", "m_BD_11"
)
jag <- rjags::jags.model(
  file = paste0("EXNEX.bugs"),
  data = data_jags,
  n.chains = n_MCMC_chain, n.adapt = n.adapt
)
update(jag, n.burnin) #may need to add burn-in phase after adaptation phase, but not always necessary
posterior_sample <- rjags::coda.samples(  #multiple chains as mcmc.list object
  model = jag,
  variable.names = pars_jags,
  n.iter = MCMC_SAMPLE,
  thin = n.thin
)
draws = as.data.frame( #single set of draws from all chains. Inference is based on all draws assuming convergence
  rbind(posterior_sample[[1]], posterior_sample[[2]], posterior_sample[[3]]))

###### RESULTS ######
setwd(outdir)
#Rhat -- have the chains converged? Rhat < 1.1 suggests convergence
R_hat = coda::gelman.diag(posterior_sample)$psrf[,1] #Rhat
#neff -- measure of autocorrelation in chains. One recommendation is neff > 100
n_eff = coda::effectiveSize(posterior_sample) #n.eff, summed across chains
Quantity = c("w_A_0", "w_A_1", "w_AC_00", "w_AC_01", "w_AC_10", "w_AC_11", "w_AD_00", "w_AD_01", "w_AD_10", "w_AD_11", 
             "w_B_0", "w_B_1", "w_BC_00", "w_BC_01", "w_BC_10", "w_BC_11", "w_BD_00", "w_BD_01", "w_BD_10", "w_BD_11", 
             "delta1", "delta2", "omega", "tau_NR", "tau_R", 
             "AAC00", "AAC01", "AAC10", "AAC11", "AAD00", "AAD01", "AAD10", "AAD11", "BBC00", "BBC01", "BBC10", "BBC11", "BBD00", "BBD01", "BBD10", "BBD11")
diag_tbl = data.frame(Quantity, R_hat, n_eff, row.names = NULL)
write.csv(diag_tbl, paste0("EXNEX_diag_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)
param_hat = apply(draws[,1:25], 2, mean)
param_var = apply(draws[,1:25], 2, var)
param_names = c("w_A_0", "w_A_1", "w_AC_00", "w_AC_01", "w_AC_10", "w_AC_11", "w_AD_00", "w_AD_01", "w_AD_10", "w_AD_11", 
                "w_B_0", "w_B_1", "w_BC_00", "w_BC_01", "w_BC_10", "w_BC_11", "w_BD_00", "w_BD_01", "w_BD_10", "w_BD_11", 
                "delta1", "delta2", "omega", "tau_NR", "tau_R")
# The parameter estimates in this case are the pathway level random effects -- not the estimates of the pathway means although those can easily be derived
param_tbl = data.frame(param_names, param_hat, sqrt(param_var), row.names = NULL)
colnames(param_tbl) = c("Parameter", "Estimate", "Standard Error")
write.csv(param_tbl, paste0("EXNEX_param_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)
DTR_hat = apply(draws[,26:41], 2, mean)
DTR_var = apply(draws[,26:41], 2, var)
DTR_quantiles = cbind(apply(draws[,26:41], 2, function(x) quantile(x, probs = 0.025)), 
                      apply(draws[,26:41], 2, function(x) quantile(x, probs = 0.975)))
DTR_names = c("AAC00", "AAC01", "AAC10", "AAC11", "AAD00", "AAD01", "AAD10", "AAD11", "BBC00", "BBC01", "BBC10", "BBC11", "BBD00", "BBD01", "BBD10", "BBD11")
DTR_tbl = data.frame(DTR_names, DTR_hat, sqrt(DTR_var), DTR_quantiles)
colnames(DTR_tbl) = c("DTR", "Estimate", "Standard Error", "95% CI - Lower Bound", "95% CI - Upper Bound")
write.csv(DTR_tbl, paste0("EXNEX_DTR_tbl_Scenario_", scenario, "_Effect_", type, "_", size, "_N_", N, ".csv"), row.names = FALSE)