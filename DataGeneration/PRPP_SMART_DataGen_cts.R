###########################################################################################################################

# Simulation code for SMART Partially Randomized Preference Trial (PRPP-SMART)

# Author: Mari Wank
# Modified by Sarah Medley for continuous outcomes
###########################################################################################################################
# Description: 
#    Code simulates a dataset from a simple 2 stage patient preference SMART design with a continuous outcome. 
#    Where applicable code is adapted from Esserman et al. 2017.
#    Please see PRPP_SMART_Data_Generation.pdf for more details on how data was simulated. 
#    Note that the below code uses a trial pathway approach for generating Y and a logit approach 
#    for modeling preference. 

###########################################################################################################################

#Function: gendata


#Purpose: This function generates one simulation dataset for a two stage PRPP-SMART trial with a binary outcome
#       The dataset has the following variables:
#       (1)ID: Numeric subject ID variable 
#       (2)X1: A continuous baseline variable generated from N(0,1)
#       (3)X2: A binary baseline variable generated from Bern(0.5). This variable represents home vs clinic/doctor’s office for treatment (from Yang et al. 2022)
#       (4)P1: Categorical variable; Preference for first stage treatment; takes values A, B, NP; A-Prefer A, B-Prefer B, NP-No preference
#       (5)S1_Preference: Binary Variable: Indicator for whether an individual had a preference in stage 1. 1-had a preference, 0-had no preference
#       (5)T1: Binary variable; Individuals treatment assigned at Stage 1; takes Values A or B. 
#       (6)R: Binary variable; Individual's stage 1 response. 1 denotes response, 0 denotes no response
#       (7)P2: Categorical variable; Preference for second stage treatment; takes values C, D, NP, 999. C-Prefer C, D-Prefer D, NP-No preference, 999-stage 1 responders are not asked stage 2 preference
#       (8)S2_Preference: Binary Variable: Indicator for whether an individual had a preference in stage 2. 1-had a preference, 0-had no preference, 999-responder in stage 1
#       (9)T2: Categorical variable; Individuals treatment assigned at Stage 2; takes values C, D, 999 (stage 1 resopnder continue on initial treatment) 
#       (10)Y: Binary variable; Individual's trial outcome. 1 denotes response, 0 denotes no response
#       (11)Trial_Path: Trial path the patient belongs to


# Required Parameters: 
#         N: Number of individuals in the generated dataset
#         pNP_target: desired proportion of individuals expressing No Preference in stage 1
#         pTheta_target: desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
#         pNP2_target: desired proportion of patients expressing No Preference in stage 2 (among non-responders)
#         pTheta2_target: desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)
#         Pa: Probability of responding given treatment A
#         Pb: Probability of responding given treatment B
#         muPAR: Probability of outcome Y for trial path 1: Preferred A in stage 1, responded
#         muRAR: Probability of outcome Y for trial path 2: Randomized A in stage 1, responded
#         muPBR: Probability of outcome Y for trial path 3: Preferred B in stage 1, responded
#         muRBR: Probability of outcome Y for trial path 4: Randomized B in stage 1, responded
#         muPANRPC: Probability of outcome Y for trial path 5: Preferred A in stage 1, no response, preferred C in stage 2
#         muPANRRC: Probability of outcome Y for trial path 6: Preferred A in stage 1, no response, randomized C in stage 2
#         muPANRPD: Probability of outcome Y for trial path 7: Preferred A in stage 1, no response, preferred D in stage 2
#         muPANRRD: Probability of outcome Y for trial path 8: Preferred A in stage 1, no response, randomized D in stage 2
#         muRANRPC: Probability of outcome Y for trial path 9: Randomized A in stage 1, no response, preferred C in stage 2
#         muRANRRC: Probability of outcome Y for trial path 10: Randomized A in stage 1, no response, randomized C in stage 2
#         muRANRPD: Probability of outcome Y for trial path 11: Randomized A in stage 1, no response, preferred D in stage 2
#         muRANRRD: Probability of outcome Y for trial path 12: Randomized A in stage 1, no response, randomized D in stage 2
#         muPBNRPC: Probability of outcome Y for trial path 13: Preferred B in stage 1, no response, preferred C in stage 2
#         muPBNRRC: Probability of outcome Y for trial path 14: Preferred B in stage 1, no response, randomized C in stage 2
#         muPBNRPD: Probability of outcome Y for trial path 15: Preferred B in stage 1, no response, preferred D in stage 2
#         muPBNRRD: Probability of outcome Y for trial path 16: Preferred B in stage 1, no response, randomized D in stage 2
#         muRBNRPC: Probability of outcome Y for trial path 9: Randomized B in stage 1, no response, preferred C in stage 2
#         muRBNRRC: Probability of outcome Y for trial path 10: Randomized B in stage 1, no response, randomized C in stage 2
#         muRBNRPD: Probability of outcome Y for trial path 11: Randomized B in stage 1, no response, preferred D in stage 2
#         muRBNRRD: Probability of outcome Y for trial path 12: Randomized B in stage 1, no response, randomized D in stage 2
#         sigma2: Data variance in each pathway. In other words, the noise in Y or variation that is not explained by the mean model.

#Output:   
#       A dataset

###########################################################################################################################

generate_data <- function(N=N, pNP_target=pNP_target, pTheta_target=pTheta_target, pNP2_target=pNP2_target, pTheta2_target=pTheta2_target, Pa=Pa, Pb=Pb, Pa1=Pa1, Pb1=Pb1, muPAR=muPAR, muRAR=muRAR, muPBR=muPBR, muRBR=muRBR, muPANRPC=muPANRPC, muPANRRC=muPANRRC, muPANRPD=muPANRPD, muPANRRD=muPANRRD, muRANRPC=muRANRPC, muRANRRC=muRANRRC, muRANRPD=muRANRPD, muRANRRD=muRANRRD, muPBNRPC=muPBNRPC, muPBNRRC=muPBNRRC, muPBNRPD=muPBNRPD,
                          muPBNRRD=muPBNRRD, muRBNRPC=muRBNRPC, muRBNRRC=muRBNRRC, muRBNRPD=muRBNRPD, muRBNRRD=muRBNRRD, sigma2=sigma2){
  
  expit <- function(y) exp(y)/(1+exp(y))
  
  
  #Function: preference 
  #Purpose: These functions assign the preference of each individual based on the true propensity for exhibiting 
  #         a preference for treatment A, B, or having no preference in stage 1 and treatment C, D, or having no
  #         no preference in stage 2. Note, pref_probx is a Nx3 matrix with the
  #         true probabilities of preference for each individual
  
  
  preference_stage1 <- function(i){
    sample(c("A", "B", "NP"), size=1, replace = TRUE, prob=pref_prob1[i,])
  }
  
  preference_stage2 <- function(i){
    sample(c("C", "D", "NP"), size=1, replace = TRUE, prob=pref_prob2[i,])
  }
  
  # load libraries 
  library(Rlab)
  library(tidyverse)
  #library(WeightIt)
  
  N <- N # number of subjects in population 
  
  # generate baseline covariates 
  X1<-rnorm(N,mean=0,sd=1)  # generate X1
  X2<-rbinom(N,1,0.5)  # generate X2: Home vs clinic/doctor’s office for treatment
  
  # set target proportion of patients expressing No Preference in stage 1
  pNP_target <- pNP_target
  # Set theta target proportion: proportion who prefer A among those with a preference (1-pNP_target)
  pTheta_target <- pTheta_target
  
  # set target proportion of patients expressing No Preference in stage 2
  pNP2_target <- pNP2_target
  # set target proportion of patients expressing preference for treatment C among those with a preference
  pTheta2_target <- pTheta2_target
  
  
  ### Preference Model: Stage 1 Propensity Scores for Preference ###
  #    Here we generate preference for first stage treatment for each patient. The true propensity for exhibiting a 
  #    preference for treatment A, B, or having no preference for each subject in the simulation population 
  #    is modeled using a logit model approach where first-stage propensity scores for preference are conditional on the 
  #    first-stage baseline covariates X1 and X2 (see PRPP_SMART_Data_Generation.pdf for details). In our preference model, 
  #    we search for an intercept via the uniroot function in r which allows us to control the proportion of subjects exhibiting 
  #    no preference/preference for A/B
  
  
  ## No Preference ##
  
  # Coefficients on X1 and X2 
  a1 <- 0.2
  a2 <- 0.133 # from Yang et al. 2022
  
  # Find intercept to get probability of no preference close to target, on average
  search0_NP <- function(a0)
  {
    alpha <- rbind(a0, a1, a2)
    lnp <- cbind(1, X1, X2) %*% alpha # log scale
    pNP <- expit(lnp) # probability scale
    mean(pNP-pNP_target) # want this to be close to 0
  }
  a0_star <- uniroot(search0_NP, c(-4, 4))$root # intercept value that will produce a marginal no preference allocation close to target
  alpha <- c(a0_star, a1, a2)
  prob_NP <- expit(cbind(1, X1, X2) %*% alpha) # calculate no preference probability for each patient
  
  
  ## Model Theta: Prefer A among those with a preference ##
  pA_marginal <- pTheta_target*(1-pNP_target) # marginal A probability in simulated data. COMPARE TO THIS 
  
  # Coefficients on X1 and X2 
  b1 <- 0.05
  b2 <- 0.1
  
  # Find intercept to get theta close to target, on average
  search0_Theta <- function(b0)
  {
    beta <- rbind(b0, b1, b2)
    lptheta<- cbind(1, X1, X2) %*% beta # log scale
    ptheta <- expit(lptheta) # probability scale
    mean(ptheta-pTheta_target) # want this to be close to 0
  }
  
  b0_star <- uniroot(search0_Theta, c(-4,4))$root 
  beta <- c(b0_star, b1, b2)
  prob_A <- expit(cbind(1, X1, X2) %*% beta) * (1-prob_NP)
  prob_B <- (1-expit(cbind(1, X1, X2) %*% beta)) * (1-prob_NP) 
  
  # prob_B <- 1-prob_NP-prob_A
  
  ## generate preference of first stage treatment of each individual ##
  
  # put probabilities into a matrix
  pref_prob1 <- cbind(prob_A,prob_B,prob_NP) 
  colnames(pref_prob1) <- c("Prob Prefer A", "Prob Prefer B", "Prob No Preference") # each row is a subject's probability to prefer A, prefer B, or have no preference
  #apply(pref_prob1,1,sum) # check
  # sample based on probabilities 
  P1 <- sapply(1:N, preference_stage1) # stage 1 preference 
  # check if allocations are close to specified targets
  #prop.table(table(P1))
  #sum(prop.table(table(P1))) # check to make sure A,B,NP sum to 1
  #length(which(P1 == "A"))/ length(which(P1 == "A" | P1 == "B")) # should be close to theta target value
  
  
  
  ### Stage 1 treatment assigned ###
  #    For each subject in the simulated population we generate their stage 1 assigned treatment
  #    based on the true propensities for exhibiting a preference for treatment A, B, or having no preference
  #    Prefer A: get A
  #    Prefer B: get B
  #    No preference: equal probability of getting treatment A/B
  
  
  # generate actual stage 1 treatment that each individual received
  # 1-treatment A, 0-treatment B
  T1 <- 
    (P1 == "A") * 1 + # prefer A get A
    (P1 == "B") * 0 + # prefer B get B
    (P1 == "NP") * sample(c(1, 0), N, replace = T, prob=c(0.5,0.5)) # no preference randomly assign A/B
  
  ### Generate Preference Indicator ###
  #    Here we create an indicator, S1_P, for whether an individual had a preference in stage 1 or not. Specifically, 1
  #    if an individual had a preference in stage 1, 0 if an individual had no preference in stage 1
  
  S1_P <- ifelse(P1 == "A" | P1 == "B", 1, 0)
  
  ### Generate Stage 1 Response variable ###
  #    For each subject we generate a binary response where 1 indicates the subjects responds
  #    to the assigned Stage 1 treatment and 0 indicates the subject does not response to the assigned 
  #    stage 1 treatment. Variable is generated by setting prespecified probabilities of responding to
  #    treatment A and B.
  
  Pi_A <- Pa # Pr(R=1|T1=A,P1=0)
  Pi_B <- Pb # Pr(R=1|T1=B,P1=0)
  Pi_A1 <- Pa1 # Pr(R=1|T1=A,P1=1)
  Pi_B1 <- Pb1 # Pr(R=1|T1=B,P1=1)
  
  # Generate stage 1 response variable 
  R<-rep(NA,N) # 0-no response to treatment, 1-responds to treatment
  
  # response paths
  R[S1_P == 1 & T1==1]<-rbinom(sum(S1_P == 1 & T1==1),1,Pi_A1) # got preferred treatment A in Stage 1
  R[S1_P == 0 & T1==1]<-rbinom(sum(S1_P == 0 & T1==1),1,Pi_A) # got randomized treatment A in Stage 1
  R[S1_P == 1 & T1==0]<-rbinom(sum(S1_P == 1 & T1==0),1,Pi_B1) # got preferred treatment B in Stage 1
  R[S1_P == 0 & T1==0]<-rbinom(sum(S1_P == 0  & T1==0),1,Pi_B) # got randomized treatment B in Stage 1
  
  ### Prefernce Model: Stage 2 Propensity Scores for Preference ###
  #    Here we generate preference of second stage treatment only for first-stage non-responders. Patients who respond
  #    continue on the stage 1 treatment. The true propensity for exhibiting a 
  #    preference for treatment C, D, or having no preference for non-responders in the simulation population 
  #    is modeled using a logit model approach where second-stage propensity scores for preference are conditional on whether
  #    the patient had a preference in stage 1  (see PRPP_SMART_Data_Generation.pdf for details). In our preference model, 
  #    we search for an intercept via the uniroot function in r which allows us to control the proportion of non-responder subjects 
  #    exhibiting no preference/preference for C/D
  
  # find non-responders 
  nr_index <- which(R == 0)
  
  ## No Preference ##
  
  # Coefficients on S1_P 
  c1 <- 0.3
  
  # Find intercept to get probability of no preference close to target, on average
  search0_NP2 <- function(c0){
    alpha <- rbind(c0, c1)
    lnp <- cbind(1, S1_P[nr_index]) %*% alpha # log scale
    pNP <- expit(lnp) # probability scale
    mean(pNP-pNP2_target) # want this to be close to 0
  }
  c0_star <- uniroot(search0_NP2, c(-4, 4))$root # intercept value that will produce a marginal no preference allocation close to target
  gamma <- c(c0_star, c1)
  prob_NP2 <- expit(cbind(1, S1_P[nr_index]) %*% gamma) # calculate no preference probability for each patient
  
  ## Model Theta: Prefer C among those with a preference ##
  pC_marginal <- pTheta2_target*(1-pNP2_target) # marginal C probability in simulated data. COMPARE TO THIS 
  
  # Coefficients on S1_P
  d1 <- 0.1
  
  # Find intercept to get theta close to target, on average
  search0_Theta2 <- function(d0)
  {
    beta <- rbind(d0, d1)
    lptheta<- cbind(1, S1_P[nr_index]) %*% beta # log scale
    ptheta <- expit(lptheta) # probability scale
    mean(ptheta-pTheta2_target) # want this to be close to 0
  }
  
  d0_star <- uniroot(search0_Theta2, c(-4,4))$root 
  phi <- c(d0_star, d1)
  prob_C <- expit(cbind(1, S1_P[nr_index]) %*% phi) * (1-prob_NP2)
  prob_D <- (1-expit(cbind(1, S1_P[nr_index]) %*% phi)) * (1-prob_NP2) 
  
  # prob_D <- 1-prob_NP2-prob_C
  
  
  ## generate preference of second stage treatment for non-responders ##
  
  # put probabilities into a matrix
  pref_prob2 <- cbind(prob_C,prob_D,prob_NP2) 
  colnames(pref_prob2) <- c("Prob Prefer C", "Prob Prefer D", "Prob No Preference") # each row is a subject's probability to prefer C, prefer D, or have no preference
  #apply(pref_prob2,1,sum) # check
  # sample based on probabilities 
  n2 <- length(nr_index)
  P2 <- sapply(1:n2, preference_stage2) # stage 2 preference only generated for stage 1 non-responders
  # check if allocations are close to specified targets
  #prop.table(table(P2))
  #sum(prop.table(table(P2))) # check to make sure A,B,NP sum to 1
  #length(which(P2 == "C"))/ length(which(P2 == "C" | P2 == "D")) # should be close to theta target value
  
  # create final stage 2 preference variable 
  P2_final <- rep(999, N) # 999 for stage 1 responders since they continue on stage 1 treatment 
  
  P2_final[nr_index] <- P2
  
  ### Stage 2 treatment assigned ###
  #    For each subject in the simulated population we generate their stage 2 assigned treatment.
  #    For non-responders this is based off the true propensities for exhibiting a preference for treatment C, D, or having no preference
  #    For responders they continue on their stage 1 assigned treatment.
  #    Non-responders:
  #      Prefer C: get C
  #      Prefer D: get D
  #      No preference: equal probability of getting treatment C/D
  
  
  # generate actual stage 2 treatment that each individual received
  # 1-treatment C, 0-treatment D
  T2 <- c() 
  for(i in 1:N){
    
    if (R[i]==1){ # responds to Stage 1 treatment
      T2[i] <- T1[i] # continues on initial treatment received 
      
    }
    if (R[i] == 0 & P2_final[i] == "NP"){# no response to Stage 1 treatment & has no preference in stage 2
      T2[i] <- rbern(1,0.5) # equal chance for treatment C/D
    }
    
    if (R[i] == 0 &  P2_final[i] == "C"){# no response to Stage 1 treatment & prefer C
      T2[i] <- 1 # get preferred treatment 
    }
    
    if (R[i] == 0 &  P2_final[i] == "D"){# no response to Stage 1 treatment & prefer D
      T2[i] <- 0 # get preferred treatment 
    }
  }
  
  
  ### Generate Stage 2 Preference Indicator ###
  #    Here we create an indicator, S2_P, for whether an individual had a preference in stage 2 or not. Specifically, 1
  #    if an individual had a preference in stage 2, 0 if an individual had no preference in stage 2
  
  S2_P <- ifelse(P2_final == "C" | P2_final == "D", 1, 0)
  
  S2_P[which(R==1)] <- 999
  
  # Relabel T1 and T2 with actual treatment 
  T1[T1==1] <- "A"
  T1[T1==0] <- "B"
  
  T2[which(R ==1)] <- 999 # assign stage 1 responders 999 
  T2[T2==1] <- "C"
  T2[T2==0] <- "D"
  
  ### Generate Trial Outcome Variable ###
  #    For each subject we generate the final binary trial outcome given their study path where 1 indicates the subjects responds
  #    0 indicates the subject does not respond.
  
  muPAR <- muPAR 
  muRAR <- muRAR
  muPBR <- muPBR
  muRBR <- muRBR
  muPANRPC <- muPANRPC
  muPANRRC <- muPANRRC
  muPANRPD <- muPANRPD
  muPANRRD <- muPANRRD
  muRANRPC <- muRANRPC
  muRANRRC <- muRANRRC
  muRANRPD <- muRANRPD
  muRANRRD <- muRANRRD
  muPBNRPC <- muPBNRPC
  muPBNRRC <- muPBNRRC
  muPBNRPD <- muPBNRPD
  muPBNRRD <- muPBNRRD
  muRBNRPC <- muRBNRPC
  muRBNRRC <- muRBNRRC
  muRBNRPD <- muRBNRPD
  muRBNRRD <- muRBNRRD
  sigma2 <- sigma2 #input sd into rnorm, not variance
  
  Y<-rep(NA,N) # 0-no response to treatment, 1-responds to treatment
  
  # response paths
  Y[S1_P == 1 & T1=="A" & R==1]<-rnorm(sum(S1_P == 1 & T1=="A" & R==1), muPAR, sqrt(sigma2))
  Y[S1_P == 0 & T1=="A" & R==1]<-rnorm(sum(S1_P == 0 & T1=="A" & R==1), muRAR, sqrt(sigma2))
  Y[S1_P == 1 & T1=="B" & R==1]<-rnorm(sum(S1_P == 1 & T1=="B" & R==1), muPBR, sqrt(sigma2))
  Y[S1_P == 0 & T1=="B" & R==1]<-rnorm(sum(S1_P == 0  & T1=="B" & R==1), muRBR, sqrt(sigma2))
  
  # paths that begin with prefer A in stage 1
  Y[S1_P == 1 & T1=="A" & R==0 & S2_P == 1 & T2 == "C"]<-rnorm(sum(S1_P == 1 & T1=="A" & R==0 & S2_P == 1 & T2 == "C"), muPANRPC, sqrt(sigma2))
  Y[S1_P == 1 & T1=="A" & R==0 & S2_P == 0 & T2 == "C"]<-rnorm(sum(S1_P == 1 & T1=="A" & R==0 & S2_P == 0 & T2 == "C"), muPANRRC, sqrt(sigma2))
  Y[S1_P == 1 & T1=="A" & R==0 & S2_P == 1 & T2 == "D"]<-rnorm(sum(S1_P == 1 & T1=="A" & R==0 & S2_P == 1 & T2 == "D"), muPANRPD, sqrt(sigma2))
  Y[S1_P == 1 & T1=="A" & R==0 & S2_P == 0 & T2 == "D"]<-rnorm(sum(S1_P == 1 & T1=="A" & R==0 & S2_P == 0 & T2 == "D"), muPANRRD, sqrt(sigma2))
  
  # paths that begin with randomize A in stage 1
  Y[S1_P == 0 & T1=="A" & R==0 & S2_P == 1 & T2 == "C"]<-rnorm(sum(S1_P == 0 & T1=="A" & R==0 & S2_P == 1 & T2 == "C"), muRANRPC, sqrt(sigma2))
  Y[S1_P == 0 & T1=="A" & R==0 & S2_P == 0 & T2 == "C"]<-rnorm(sum(S1_P == 0 & T1=="A" & R==0 & S2_P == 0 & T2 == "C"), muRANRRC, sqrt(sigma2))
  Y[S1_P == 0 & T1=="A" & R==0 & S2_P == 1 & T2 == "D"]<-rnorm(sum(S1_P == 0 & T1=="A" & R==0 & S2_P == 1 & T2 == "D"), muRANRPD, sqrt(sigma2))
  Y[S1_P == 0 & T1=="A" & R==0 & S2_P == 0 & T2 == "D"]<-rnorm(sum(S1_P == 0 & T1=="A" & R==0 & S2_P == 0 & T2 == "D"), muRANRRD, sqrt(sigma2))
  
  # paths that begin with prefer B in stage 1
  Y[S1_P == 1 & T1=="B" & R==0 & S2_P == 1 & T2 == "C"]<-rnorm(sum(S1_P == 1 & T1=="B" & R==0 & S2_P == 1 & T2 == "C"), muPBNRPC, sqrt(sigma2))
  Y[S1_P == 1 & T1=="B" & R==0 & S2_P == 0 & T2 == "C"]<-rnorm(sum(S1_P == 1 & T1=="B" & R==0 & S2_P == 0 & T2 == "C"), muPBNRRC, sqrt(sigma2))
  Y[S1_P == 1 & T1=="B" & R==0 & S2_P == 1 & T2 == "D"]<-rnorm(sum(S1_P == 1 & T1=="B" & R==0 & S2_P == 1 & T2 == "D"), muPBNRPD, sqrt(sigma2))
  Y[S1_P == 1 & T1=="B" & R==0 & S2_P == 0 & T2 == "D"]<-rnorm(sum(S1_P == 1 & T1=="B" & R==0 & S2_P == 0 & T2 == "D"), muPBNRRD, sqrt(sigma2))
  
  
  # paths that begin with randomize B in stage 1
  Y[S1_P == 0 & T1=="B" & R==0 & S2_P == 1 & T2 == "C"]<-rnorm(sum(S1_P == 0 & T1=="B" & R==0 & S2_P == 1 & T2 == "C"), muRBNRPC, sqrt(sigma2))
  Y[S1_P == 0 & T1=="B" & R==0 & S2_P == 0 & T2 == "C"]<-rnorm(sum(S1_P == 0 & T1=="B" & R==0 & S2_P == 0 & T2 == "C"), muRBNRRC, sqrt(sigma2))
  Y[S1_P == 0 & T1=="B" & R==0 & S2_P == 1 & T2 == "D"]<-rnorm(sum(S1_P == 0 & T1=="B" & R==0 & S2_P == 1 & T2 == "D"), muRBNRPD, sqrt(sigma2))
  Y[S1_P == 0 & T1=="B" & R==0 & S2_P == 0 & T2 == "D"]<-rnorm(sum(S1_P == 0 & T1=="B" & R==0 & S2_P == 0 & T2 == "D"), muRBNRRD, sqrt(sigma2))
  
  ### Create final dataframe ### 
  
  
  ### Create trial path variable ###
  trial_path <- rep(NA, N)
  
  # responder paths
  trial_path[which(S1_P==1 & T1=="A" & R==1)] <- "PAR"
  trial_path[which(S1_P==0 & T1=="A" & R==1)] <- "RAR"
  trial_path[which(S1_P==1 & T1=="B" & R==1)] <- "PBR"
  trial_path[which(S1_P==0 & T1=="B" & R==1)] <- "RBR"
  
  # non-responder paths
  trial_path[which(S1_P==1 & T1=="A" & R==0 & S2_P==1 & T2=="C")] <- "PANRPC"
  trial_path[which(S1_P==1 & T1=="A" & R==0 & S2_P==0 & T2=="C")] <- "PANRRC"
  trial_path[which(S1_P==1 & T1=="A" & R==0 & S2_P==1 & T2=="D")] <- "PANRPD"
  trial_path[which(S1_P==1 & T1=="A" & R==0 & S2_P==0 & T2=="D")] <- "PANRRD"
  
  trial_path[which(S1_P==1 & T1=="B" & R==0 & S2_P==1 & T2=="C")] <- "PBNRPC"
  trial_path[which(S1_P==1 & T1=="B" & R==0 & S2_P==0 & T2=="C")] <- "PBNRRC"
  trial_path[which(S1_P==1 & T1=="B" & R==0 & S2_P==1 & T2=="D")] <- "PBNRPD"
  trial_path[which(S1_P==1 & T1=="B" & R==0 & S2_P==0 & T2=="D")] <- "PBNRRD"
  
  trial_path[which(S1_P==0 & T1=="A" & R==0 & S2_P==1 & T2=="C")] <- "RANRPC"
  trial_path[which(S1_P==0 & T1=="A" & R==0 & S2_P==0 & T2=="C")] <- "RANRRC"
  trial_path[which(S1_P==0 & T1=="A" & R==0 & S2_P==1 & T2=="D")] <- "RANRPD"
  trial_path[which(S1_P==0 & T1=="A" & R==0 & S2_P==0 & T2=="D")] <- "RANRRD"
  
  trial_path[which(S1_P==0 & T1=="B" & R==0 & S2_P==1 & T2=="C")] <- "RBNRPC"
  trial_path[which(S1_P==0 & T1=="B" & R==0 & S2_P==0 & T2=="C")] <- "RBNRRC"
  trial_path[which(S1_P==0 & T1=="B" & R==0 & S2_P==1 & T2=="D")] <- "RBNRPD"
  trial_path[which(S1_P==0 & T1=="B" & R==0 & S2_P==0 & T2=="D")] <- "RBNRRD"
  
  data_output <- 
    tibble(
      id = 1:N,
      X1 = X1,
      X2 = X2,
      P1 = P1,
      S1_Preference = S1_P,
      T1 = T1,
      R = R,
      P2 = P2_final,
      S2_Preference = S2_P,
      T2 = T2,
      Y = Y,
      Trial_Path = trial_path,
    )
  
  # create list of sample sizes for each DTR
  DTR_n = c()
  DTR_n[1] = sum(data_output$Trial_Path == "RAR") + 
    sum(data_output$Trial_Path == "RANRRC") #AC00
  DTR_n[2] = sum(data_output$Trial_Path == "RAR") +
    sum(data_output$Trial_Path == "RANRRD") #AD00
  DTR_n[3] = sum(data_output$Trial_Path == "RBR") + 
    sum(data_output$Trial_Path == "RBNRRC") #BC00
  DTR_n[4] = sum(data_output$Trial_Path == "RBR") +
    sum(data_output$Trial_Path == "RBNRRD") #BD00
  DTR_n[5] = sum(data_output$Trial_Path == "RAR") + 
    sum(data_output$Trial_Path == "RANRPC") #AC01
  DTR_n[6] = sum(data_output$Trial_Path == "RAR") +
    sum(data_output$Trial_Path == "RANRPD") #AD01
  DTR_n[7] = sum(data_output$Trial_Path == "RBR") + 
    sum(data_output$Trial_Path == "RBNRPC") #BC01
  DTR_n[8] = sum(data_output$Trial_Path == "RBR") +
    sum(data_output$Trial_Path == "RBNRPD") #BD01
  DTR_n[9] = sum(data_output$Trial_Path == "PAR") + 
    sum(data_output$Trial_Path == "PANRRC") #AC10
  DTR_n[10] = sum(data_output$Trial_Path == "PAR") +
    sum(data_output$Trial_Path == "PANRRD") #AD10
  DTR_n[11] = sum(data_output$Trial_Path == "PBR") + 
    sum(data_output$Trial_Path == "PBNRRC") #BC10
  DTR_n[12] = sum(data_output$Trial_Path == "PBR") +
    sum(data_output$Trial_Path == "PBNRRD") #BD10
  DTR_n[13] = sum(data_output$Trial_Path == "PAR") + 
    sum(data_output$Trial_Path == "PANRPC") #AC11
  DTR_n[14] = sum(data_output$Trial_Path == "PAR") +
    sum(data_output$Trial_Path == "PANRPD") #AD11
  DTR_n[15] = sum(data_output$Trial_Path == "PBR") + 
    sum(data_output$Trial_Path == "PBNRPC") #BC11
  DTR_n[16] = sum(data_output$Trial_Path == "PBR") +
    sum(data_output$Trial_Path == "PBNRPD") #BD11
  
  
  # relabel data 
  data_output_relabel <- data_output %>% mutate(T1 = replace(T1, T1 == "A", 1),
                                                T1 = replace(T1, T1 == "B", 0),
                                                T1 = as.numeric(T1),
                                                T2 = replace(T2, T2 == "C", 1),
                                                T2 = replace(T2, T2 == "D", 0),
                                                T2 = replace(T2, T2 == "999", 0),
                                                T2 = as.numeric(T2),
                                                S2_Preference = replace(S2_Preference, S2_Preference == 999, 0),
                                                S2_Preference = as.numeric(S2_Preference))
  
  
  ## Replicate and weight Data - need four datasets because responders consistent with four DTRs
  # get empirical proportions for weight calculation
  pNP_target_hat <- sum(data_output_relabel$S1_Preference == 0)/nrow(data_output_relabel)
  pNP2_target_hat <- sum(data_output_relabel$S2_Preference == 0 & data_output_relabel$R == 0)/sum(data_output_relabel$R == 0)
  pNP2_P1_1 <- sum(data_output_relabel$S2_Preference == 0 & data_output_relabel$R == 0 & data_output_relabel$S1_Preference == 1)/sum(data_output_relabel$R == 0 & data_output_relabel$S1_Preference == 1)
  pNP2_P1_0 <- sum(data_output_relabel$S2_Preference == 0 & data_output_relabel$R == 0 & data_output_relabel$S1_Preference == 0)/sum(data_output_relabel$R == 0 & data_output_relabel$S1_Preference == 0)
  pTheta_target_hat <- sum(data_output_relabel$T1 == 1 & data_output_relabel$S1_Preference == 1) / sum(data_output_relabel$S1_Preference == 1)
  pTheta2_target_hat <- sum(data_output_relabel$T2 == 1 & data_output_relabel$R == 0 & data_output_relabel$S2_Preference == 1) / sum(data_output_relabel$S2_Preference == 1)
  pTheta2_P1_1 <- sum(data_output_relabel$T2 == 1 & data_output_relabel$R == 0 & data_output_relabel$S2_Preference == 1 & data_output_relabel$S1_Preference == 1) / sum(data_output_relabel$S2_Preference == 1 & data_output_relabel$S1_Preference == 1)
  pTheta2_P1_0 <- sum(data_output_relabel$T2 == 1 & data_output_relabel$R == 0 & data_output_relabel$S2_Preference == 1 & data_output_relabel$S1_Preference == 0) / sum(data_output_relabel$S2_Preference == 1 & data_output_relabel$S1_Preference == 0)
  
  # Inverse Probability of Treatment Weights 
  # P(P1)*P(T1 | P1)*P(P2 | T1, P1)*P(T2 | P1, T1, P2)
  # Since the generation of P2 and T2 only depends on P1, this simplifies here to
  # P(P1)*P(T1 | P1)*P(P2 | P1)*P(T2 | P1)
  
  data_output_relabel$w=rep(0)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RAR")] <- 1/(pNP_target_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RANRPC")] <- 1/(pNP_target_hat*0.5*(1-pNP2_P1_0)*pTheta2_P1_0)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RANRPD")] <- 1/(pNP_target_hat*0.5*(1-pNP2_P1_0)*(1-pTheta2_P1_0))
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RANRRC")] <- 1/(pNP_target_hat*0.5*pNP2_P1_0*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RANRRD")] <- 1/(pNP_target_hat*0.5*pNP2_P1_0*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RBR")] <- 1/(pNP_target_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RBNRPC")] <- 1/(pNP_target_hat*0.5*(1-pNP2_P1_0)*pTheta2_P1_0)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RBNRPD")] <- 1/(pNP_target_hat*0.5*(1-pNP2_P1_0)*(1-pTheta2_P1_0))
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RBNRRC")] <- 1/(pNP_target_hat*0.5*pNP2_P1_0*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="RBNRRD")] <- 1/(pNP_target_hat*0.5*pNP2_P1_0*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PBR")] <- 1/(pNP_target_hat*(1-pTheta_target_hat))
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PBNRPC")] <- 1/((1-pNP_target_hat)*(1-pTheta_target_hat)*(1-pNP2_P1_1)*pTheta2_P1_1)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PBNRPD")] <- 1/((1-pNP_target_hat)*(1-pTheta_target_hat)*(1-pNP2_P1_1)*(1-pTheta2_P1_1))
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PBNRRC")] <- 1/((1-pNP_target_hat)*(1-pTheta_target_hat)*pNP2_P1_1*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PBNRRD")] <- 1/((1-pNP_target_hat)*(1-pTheta_target_hat)*pNP2_P1_1*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PAR")] <- 1/((1-pNP_target_hat)*pTheta_target_hat)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PANRPC")] <- 1/((1-pNP_target_hat)*pTheta_target_hat*(1-pNP2_P1_1)*pTheta2_P1_1)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PANRPD")] <- 1/((1-pNP_target_hat)*pTheta_target_hat*(1-pNP2_P1_1)*(1-pTheta2_P1_1))
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PANRRC")] <-1/((1-pNP_target_hat)*pTheta_target_hat*pNP2_P1_1*0.5)
  data_output_relabel$w[which(data_output_relabel$Trial_Path=="PANRRD")] <- 1/((1-pNP_target_hat)*pTheta_target_hat*pNP2_P1_1*0.5)
  
  # 1st dataset of responders setting T2=1, 
  datareps1 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps1$T2 <- 1
  datareps1$S2_Preference <- 0
  
  # 2nd dataset of responders setting T2=1, 
  datareps2 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps2$T2 <- 1
  datareps2$S2_Preference <- 1
  
  # 3rd dataset of responders setting T2=-1, 
  datareps3 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps3$T2 <- 0
  datareps3$S2_Preference <- 0
  
  # 4th dataset of responders setting T2=-1, 
  datareps4 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps4$T2 <- 0
  datareps4$S2_Preference <- 1
  
  # dataset for non-responders
  datanoresp <- data_output_relabel[data_output_relabel$R==0,]
  
  # replicated data
  replicated_dat <- rbind(datareps1,datareps2,datareps3,datareps4,datanoresp)
  
  # create data used in bayesian analysis  - change replicated_dat to replicatedweighted_dat if using upsampled data and not weights argument 
  return(list(data_output_relabel, replicated_dat, DTR_n))
}
