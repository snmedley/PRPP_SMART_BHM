## Runs more efficient version of the BHM BUGS code

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
type = 1 # 1, 2, or 3
size = "small" # "small", "moderate", or "large"
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
true_mean = expected_pref[c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)]
#Note: efficient version considers DTRs in a different order

##### Choose Prior Settings and MCMC Settings #####
setwd("../JAGSFiles")
load("Mu_Prior1.RData")
tau_sd = 3/2 # 3/2 for a prior assuming substantial heterogeneity, 3 for large heterogeneity
## JAGS settings
pars_full = c("delta1", "delta2", "dtr", "omega", "p_R", "tau_NR", "tau_R")
n.adapt <- 5000 
n.burnin <- 5000
MCMC_SAMPLE <- 10000 #number of posterior samples AFTER adaptation and burn-in
n.thin <- 1 #chains should be 2000 draws long each
n_MCMC_chain <- 3
#https://stackoverflow.com/questions/38701100/how-to-interpret-some-syntax-n-adapt-update-in-jags

## Path df and DTR df
setwd("..")
dtr_df <- read.csv(paste0("PRPP-SMART_Design_EmbeddedDTRs.csv"), header = TRUE)
n_DTRs = nrow(dtr_df)
dtr.names = dtr_df$Consistent_DTR
dtr_df$T1_idx = case_when(dtr_df$T1 == "A" & dtr_df$P1 == 0 ~ 1,
                          dtr_df$T1 == "A" & dtr_df$P1 == 1 ~ 2,
                          dtr_df$T1 == "B" & dtr_df$P1 == 0 ~ 3,
                          dtr_df$T1 == "B" & dtr_df$P1 == 1 ~ 4)
path_df <- read.csv(paste0("PRPP-SMART_Design_Pathways.csv"), header = TRUE) 
num_paths = nrow(path_df)

###### SIMULATION ######
seed = 1234 
set.seed(seed)
data <- generate_data(N=N, pNP_target=pNP_target, pTheta_target=pTheta_target, pNP2_target=pNP2_target, pTheta2_target=pTheta2_target, Pa=Pa, Pb=Pb, Pa1=Pa1, Pb1=Pb1, muPAR=muPAR, muRAR=muRAR, muPBR=muPBR, muRBR=muRBR, muPANRPC=muPANRPC, muPANRRC=muPANRRC, muPANRPD=muPANRPD, muPANRRD=muPANRRD, muRANRPC=muRANRPC, muRANRRC=muRANRRC, muRANRPD=muRANRPD, muRANRRD=muRANRRD, muPBNRPC=muPBNRPC, muPBNRRC=muPBNRRC, muPBNRPD=muPBNRPD, muPBNRRD=muPBNRRD, muRBNRPC=muRBNRPC, muRBNRRC=muRBNRRC, muRBNRPD=muRBNRPD, muRBNRRD=muRBNRRD, sigma2=sigma2)
trialpath_df <- data[[1]] %>% dplyr::group_by(Trial_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
skip = ifelse(nrow(trialpath_df) < 20, 1, 0) # skip simulation if there are positivity issues
n.DTR = data[[3]] # number of participants consistent with each DTR
df <- data[[1]] # data for prpp_smart full analysis (all subjects)
## Estimate variability in the outcome in each trial pathway
df <- df %>% group_by(Trial_Path) %>% mutate(s = sd(Y)) %>% ungroup()
df$Pathway = ifelse(df$R == 1, paste0(ifelse(df$T1 == 1, "A", "B"), "_", df$S1_Preference), 
                    paste0(ifelse(df$T1 == 1, "A", "B"), "_", df$S1_Preference, " ", ifelse(df$T2 == 1, "C", "D"), "_", df$S2_Preference))
df = left_join(df, path_df[,c(8,9)], by = "Pathway", relationship = "many-to-one") #Pathway & Path_Index columns
# Create data list for JAGS
df$g_idx = case_when(df$T1 == "1" & df$R == 1 ~ 1,
                     df$T1 == "0" & df$R == 1 ~ 2,
                     df$T1 == "1" & df$R == 0 & df$T2 == "1" ~ 3,
                     df$T1 == "1" & df$R == 0 & df$T2 == "0" ~ 4,
                     df$T1 == "0" & df$R == 0 & df$T2 == "1" ~ 5,
                     df$T1 == "0" & df$R == 0 & df$T2 == "0" ~ 6)
  
data_full = c(list(
    # PRIOR SETTINGS
    c(mu_A_p, mu_B_p, mu_AC_p, mu_AD_p, mu_BC_p, mu_BD_p), omega_p, delta1_p, delta2_p, epsilon_p, tau_sd,
    # RESPONSE DATA
    c(sum(df$T1 == 1 & df$S1_Preference == 0 & df$R == 1), sum(df$T1 == 1 & df$S1_Preference == 1 & df$R == 1), 
      sum(df$T1 == 0 & df$S1_Preference == 0 & df$R == 1), sum(df$T1 == 0 & df$S1_Preference == 1 & df$R == 1)),
    c(sum(df$T1 == 1 & df$S1_Preference == 0), sum(df$T1 == 1 & df$S1_Preference == 1), 
      sum(df$T1 == 0 & df$S1_Preference == 0), sum(df$T1 == 0 & df$S1_Preference == 1)),
    # INDIVIDUAL LEVEL DATA
    nrow(df), df$R, df$Y, df$s, df$S1_Preference, df$S2_Preference, df$Path_Index,
    # DTR LEVEL DATA
    nrow(dtr_df), dtr_df$T1_idx, dtr_df$Responder_Index, dtr_df$Nonresponder_Index,
    dtr_df$P1, dtr_df$P2,
    # TREATMENT SEQUENCE/PATHWAY LEVEL DATA
    nrow(path_df),
    length(unique(path_df$Sequence)),
    path_df$g_idx,
    (path_df$R == 1)*1
))
names(data_full) = c("mu_p", "omega_p", "delta1_p", "delta2_p", "epsilon_p", "tau_sd",
                     "N_R", "N_T1_P1",
                     "N", "R", "Y", "s","P1", "P2", "path_idx",
                     "N_dtr", "T1_idx", "R_path_idx", "NR_path_idx", "dtr_P1", "dtr_P2",
                     "N_paths", "N_trt_seq", "g_idx", "I_R")
jag <- rjags::jags.model(
    file = paste0("BHM_Efficient.bugs"),
    data = data_full,
    n.chains = n_MCMC_chain, n.adapt = n.adapt
)
# Run JAGS model -- burn-in then posterior sampling
update(jag, n.burnin) #may need to add burn-in phase after adaptation phase, but not always necessary
posterior_full <- rjags::coda.samples(  #multiple chains as mcmc.list object
    model = jag,
    variable.names = pars_full,
    n.iter = MCMC_SAMPLE,
    thin = n.thin
)
draws_full = as.data.frame( #single set of draws from all chains. Inference is based on all draws assuming convergence
    rbind(posterior_full[[1]], posterior_full[[2]], posterior_full[[3]]))
R_hat = coda::gelman.diag(posterior_full)$psrf[,1]
n_eff = coda::effectiveSize(posterior_full)
param_hat = apply(draws_full[,c(1:2,19:25)], 2, mean)
param_var = apply(draws_full[,c(1:2,19:25)], 2, var)
DTR_hat = apply(draws_full[,3:18], 2, mean)
DTR_var = apply(draws_full[,3:18], 2, var)
dtr_quantiles = cbind(apply(draws_full[,3:18], 2, function(x) quantile(x, probs = 0.025)), 
                      apply(draws_full[,3:18], 2, function(x) quantile(x, probs = 0.975)))
# Since JAGS puts parameters in alphabetical order, need to reorder true DTR means
ci_length <- dtr_quantiles[,2] - dtr_quantiles[,1]
dtr_in_ci <- data.table::between(true_mean, dtr_quantiles[,1], dtr_quantiles[,2])