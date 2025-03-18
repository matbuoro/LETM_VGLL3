rm(list=ls())   # Clear memory

# DIRECTORY ####
#setwd("~/Documents/RESEARCH/PROJECTS/VGLL3")

# PACKAGES ####
require(nimble)
library(parallel)
#library(igraph)
library(MCMCvis)
library(writexl)

# DATA ####
source("code/DATA_VGLL3_Scorff.R")
#attach(data)
#str(data)
#data$X.scaled <- NULL
dataToNimble <- list(Y=data$Y
                     #, X=data$X
                     , X.scaled=data$X.scaled
                     )
constants <- list(N=data$N
                  #, X= data$X
                  #, mean.X=data$mean.X
                  #, sd.X=data$sd.X
                  , sex=data$sex
                  #, n=data$n
                  , g=data$g
                  #, year=data$year
                  , freq=data$freq
)


# MODEL ####
#source("code/MODEL_LETM_VGLL3.R")
LETM <- nimbleCode({
    # Y: Maturation state; 1 if mature (1SW), 0 otherwise (MSW)
    # X: observable cue (growth)
    # eta: individual Proximate cue
    # theta: individual thresholds
    # a: genotypic values
    # k: scaled dominance deviation
    # sex: Male 1 / female 2
    
  #   for (s in 1:2){
  #     for (gen in 1:3){
  #       mu_X[gen,s] ~ dnorm(0,0.001)
  #       
  #       sigma_X[gen,s]~dgamma(2,1/sig)#dunif(0,10)
  #       sigma2_X[gen,s] <- pow(sigma_X[gen,s],2)
  #     }}
  # sig~dchisqr(2)
    
    #======== LIKELIHOOD ========#
    for (i in 1:N) {
      
      ## Environmental cue (X):
      #X[i]~dnorm(mu_X[g[i],sex[i]],  1/sigma2_X[g[i],sex[i]]) # Growth distributions
      #X.scaled[i] <- (X[i]-mu_X[g[i],sex[i]])/sigma_X[g[i],sex[i]] # scaled growth
      
      ## Phenotypes:
      Y[i]~dbern(p_mat[i])
      p_mat[i]<-phi(z[i]) # probit link function
      #z[i]<-((X.scaled[i] - theta[i]) / sqrt(sigma2_eta[sex[i]]))
      z[i]<-(((X.scaled[i]) / sqrt(sigma2_eta[g[i],sex[i]] + 1)) - theta[i]) / sqrt(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]] + 1))
      
      ## Threshold
      theta[i]~dnorm(mu[i,2], 1/sigma2_res[g[i],sex[i]])
      mu[i,2] <- mu_theta[sex[i]] + alpha[g[i],sex[i]]
      #mu[i,2] <- mu_theta[period[i],sex[i]] + alpha[g[i],sex[i]]
      
      ## Proximate cue (eta)
      #lower[i] <- ifelse(Y[i]==1,theta[i],-100)
      lower[i] <- Y[i] * theta[i] + (1 - Y[i]) * (-100)
      #upper[i] <- ifelse(Y[i]==0,theta[i], 100)
      upper[i] <- (1 - Y[i]) * theta[i] + Y[i] * (100)
      mu_eta[i] <- X.scaled[i] / sqrt(sigma2_eta[g[i],sex[i]] + 1) # normalized
      #eta[i] ~dnorm(mu_eta[i], 1/(sigma2_eta[g[i],sex[i]]/ (sigma2_eta[g[i],sex[i]]+ 1)));T(lower[i],upper[i]) # /!\ using model as R function, the truncated normal is not a valid R expression. You can fool the R interpreter by inserting a semicolon
      eta[i] ~T(dnorm(mu_eta[i], 1/(sigma2_eta[g[i],sex[i]]/ (sigma2_eta[g[i],sex[i]]+ 1))), lower[i],upper[i])
      
    } # End of loop i
    
    # Average thresholds
    mu_theta[1]~dnorm(0, 1) # male
    #mu_theta[2]~dnorm(0, 1) # female
    mu_theta[2] <- mu_theta[1] + delta_theta
    delta_theta~dnorm(0, 1)

    for (s in 1:2){    # sex
      # Genotype x Sex effects
      alpha[1,s] <- -a[s] # EE
      alpha[2,s] <- a[s]*k[s] # EL
      alpha[3,s] <- a[s] # LL
      
      # genotypic values
      a[s]~dnorm(0, 1)
      
      # dominance deviation
      d[s] <- k[s]*a[s] # dominance 
      k[s]~dnorm(0, 1) # scaled dominance
      
    } # end loop s
    
    
    delta_a <- a[1] - a[2]
    delta_k <- k[1] - k[2] 
    
    
    # residual variance of the latent threshold
    for (s in 1:2){
      
     mu_alpha[s] <- (a[s]*(-freq[1,s]+k[s]*freq[2,s]+freq[3,s]))
      
      for (j in 1:3){
        
        # Genotype (VGLL3) x Sex variance
        #sigma2_alpha[j,s] <- pow(alpha[j,s] - mu_alpha[s],2)*(n[j,s]/sum(n[1:3,s])) # variance genotype
        sigma2_alpha[j,s] <- pow(alpha[j,s] - mu_alpha[s],2)*(freq[j,s]) # variance genotype

      } # end loop j
    } # end loop s
    #nu <- 2
   
    for (j in 1:3){
      
      # Variance proximate cue
      #sigma2_eta[j,1] <- sigma2_X[j,1] *(ratio[1]/(1-ratio[1]))
      sigma2_eta[j,1] <- 1*(ratio[1]/(1-ratio[1])) # sigma2_Xscaled = 1
      sigma2_eta[j,2] <- 1*(ratio[2]/(1-ratio[2])) # sigma2_Xscaled = 1
      #sigma2_eta[j,2] <- sigma2_eta[j,1]
     
      # Hypothesis: same  residual  varoance by genotype and sex
      sigma2_res[j,1] <- sum(sigma2_alpha[1:3,1])*(ratio[3]/(1-ratio[3]))
      sigma2_res[j,2] <- sum(sigma2_alpha[1:3,2])*(ratio[4]/(1-ratio[4]))
      #sigma2_res[j,2] <- sigma2_res[j,1] # same for sex
      
    }
    
    # Contribution of residual variances
    # ratio[1]~dbeta(3,3)
    # ratio[2]~dbeta(3,3)
    # ratio[3]~dbeta(3,3)
    logit(ratio[1])~dnorm(0,1.5)
    logit(ratio[2])~dnorm(0,1.5)
    logit(ratio[3])~dnorm(0,1.5)
    ratio[4] <- ratio[3]
    #ratio[4]~dbeta(3,3)
    
    varE[1] <- 1# sigma2_eta[1] + 1#sum(sigma_X[1:3, 1])*(n[j,1]/sum(n[1:3,1])) # Variance environnemental

    sigma2_T[1] <- sigma2_G[1] + varE[1]
    sigma2_G[1] <-  sum(sigma2_alpha[1:3,1]) + sigma2_res[1,1] #sum(sigma2_THETA[1:3,1]) # variance genetic
    h2[1] <- sigma2_G[1] / sigma2_T[1]
    
    # female
    varE[2] <- 1#sigma2_eta[2] + 1
    sigma2_T[2] <- sigma2_G[2] + varE[2]
    sigma2_G[2] <- sum(sigma2_alpha[1:3,2]) + sigma2_res[1,2] #sum(sigma2_THETA[1:3,2]) # variance genetic
    h2[2] <- sigma2_G[2] / sigma2_T[2]
    
    
    
    
    for (s in 1:2){    # sex
      p[s] <- sum((2*freq[1,s])+(1*freq[1,s])+(0*freq[1,s]))/(2*sum(freq[1:3,s])) # population frequence allelic for E
      q[s] <- 1 - p[s]
      
      #Average effect of the gene substitution
      # Allele substitution effect of vgll3
      gamma[s] <- a[s] + d[s]*(q[s]-p[s])
      #alpha1= q[s]*gamma[s]
      #alpha2= -p[s]*gamma[s]
      
      # Heterozygosity at locus vgl3 in the case of Hardy-Weinberg equilibrium
      h[s] <- 2*p[s]*q[s]
      He[s] <- 1-(pow(p[s],2)+ pow(q[s],2))
      
      #Population mean (for VGLL3)
      M[s] <- a[s]*(p[s]-q[s])+ (2*d[s]*p[s]*q[s])
      
      # Phenotypic variance
      # varP = varG + varE = varA + varD + varI + varE
      varT[s] <- varG[s] + varE[s]
      #genetic mean of the population
      varG[s] <- varA[s] + varD[s] + sigma2_res[1,s]
      # additive genetic variance
      varA[s] <- 2*p[s]*q[s]*pow(a[s]+d[s]*(q[s]-p[s]),2)
      #varA[s] <- 2*p*q*pow(a[s]+d[s]*(q-p),2)
      # dominance genetic variance
      varD[s] <- pow(2*p[s]*q[s]*d[s],2)
      #varD[s] <- pow(2*p*q*d[s],2)
      
      # Contributions to total variance
      h2[s+2] <- varG[s] / varT[s] # total genetic variance contribution = broad-sense heritability
      h2[s+4] <- varA[s] / varT[s] # addtive genetic variance contribution = narrow-sense heritability
      h2[s+6] <- varD[s] / varG[s] # non-additive genetic variance contribution (dominance + allelic interactions)
      h2[s+8] <- varD[s] / varT[s] # non-additive genetic variance contribution (dominance + allelic interactions)
      
    }

})


# Parameters ####
parameters <- c(
  #"Y"
  #"mu_X","sigma2_X","sigma_X"
  #,"sig"
  #,"beta","eps_X"
  #,"mean_theta","mean_eta","delta_eta"
  "delta_theta"
  ,"mu_theta", "mu_alpha"#,"eps_theta"
  #,"delta_X"
  ,"delta_a"
  ,"delta_k"
  #,"delta_res"
  ,"alpha", "d","a","k"#,"k.prior"
  #,"chSq","chSq.prior"
  ,"sigma2_eta"
  #,"sigma2_e"
  ,"ratio"
  #, "sigma2_theta"
  ,"sigma2_alpha","sigma2_res"
  #,"sigma2_THETA"
  ,"h2","sigma2_T","sigma2_G"
  ,"theta","eta"
  ,"gamma","M","h","He"
  #,"gap"
  #,"delta","delta_theta","delta_alpha"
  #,"mu_alpha"
  #,"var_alpha","var_theta"
  #,"var_alpha_sex","var_alpha_gen"
  #,"sigma2_GENOTYPE", "sigma2_VGLL3","sigma2_SEX","sigma2_GENE","sigma2_TOT","sigma2_RES","h"
  ,"varG","varA","varD","varT","varE"
  #,"delta_sigma"
  # ,"p"
  
  #,"h"
  #,"psi","lpsi","psi.prior","lpsi.prior"
  #,"theta.pred"
  
  #,"delta_X_sex","delta_sdX_sex","delta_X","delta_sdX", "P_delta","smd_X"
)


# Initial values ####
mu_theta=c(0,3.8)
eta_inits <- theta_inits <- numeric(constants$N)
for (i in 1:constants$N){
  eta_inits[i] <- rnorm(1,0,1)
  if(data$Y[i]==1) {
    theta_inits[i] <- mu_theta[data$sex[i]]
    eta_inits[i] <- theta_inits[i]+0.1
  }
  if(data$Y[i]==0) {
    theta_inits[i] <- mu_theta[data$sex[i]]
    eta_inits[i] <- theta_inits[i]-0.1
  }
}

inits <- function(){
  list(
    eta=eta_inits, theta=theta_inits
    , mu_theta=c(0,NA),delta_theta=3.8
    , a=c(3,1)#,delta_a=0
    , k=c(-1,0)#,delta_k=1
    , ratio=c(rep(0.5, 3),NA)
    #, mu_X=data$mean.X
    #, sigma_X=data$sd.X
    #, sig=30
  )
}


# RUN MCMC ####
n_chains <- 2 # number of chains
n_store <- 5000 # target of number of iteration to store per chain
n_burnin <- 1000 # number of iterations to discard
n_thin <- 100 # thinning interval
n_iter <- (n_store * n_thin) + n_burnin # number of iterations to run per chain
print(n_iter)


start_time <- Sys.time()

  samples <- nimbleMCMC(code = LETM,     # model code
                            data = dataToNimble,                  # data
                            constants =constants,        # constants
                            inits = inits,          # initial values
                            monitors = parameters,   # parameters to monitor
                            WAIC=FALSE,                      #waic
                            niter = n_iter,                  # nb iterations
                            nburnin = n_burnin,              # length of the burn-in
                            nchains = n_chains,              # nb of chains
                            thin = n_thin,                   # thinning interval (default = 1)
                            samplesAsCodaMCMC=T
  )             #coda

  # Define a function to run one chain
  # run_chain <- function(seed) {
  #   nimbleMCMC(
  #     code = LETM, 
  #     data = dataToNimble, 
  #     constants = constants, 
  #     inits = inits,  # Each chain should have different inits
  #     monitors = parameters, 
  #     WAIC = FALSE, 
  #     niter = n_iter, 
  #     nburnin = n_burnin, 
  #     nchains = 1,  # Only run one chain per function call
  #     thin = n_thin, 
  #     samplesAsCodaMCMC = TRUE
  #   )
  # }
  
  # Run chains in parallel
  #samples_list <- mclapply(1:n_chains, run_chain, mc.cores = n_chains)
  # Combine chains
  #samples <- mcmc.list(samples_list)
  
# End time measurement
end_time <- Sys.time()
# Calculate and print the time difference
time_taken <- end_time - start_time
print(time_taken)

save(data,samples,time_taken, file="results/RESULTS_vgll3_scorff_nimble.RData")

#write.csv2(samples$BUGSoutput$summary, file="results/Summary_2025.csv")


# RESULTS ####

## load package
#library(MCMCvis)

# samples$chain1 <- samples$chain1[,-c(43,44)]
# samples$chain2 <- samples$chain2[,-c(43,44)]
# samples$chain3 <- samples$chain3[,-c(43,44)]

## Numerical summaries
#MCMCsummary(object = samples, round = 2)

## Trace and posterior density ####
MCMCtrace(object = samples,
          filename="MCMC_nimble_2025.pdf",
          wd = "results/",
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE, # add Rhat
          n.eff = TRUE, # add eff sample size
          params = parameters,
          excl=c("theta","eta")
)


parToPlot <-  c(
#"mu_X","sigma2_X","sigma_X"
#"sig"
"mu_theta", "mu_alpha"
,"a","k"
,"sigma2_eta"
,"ratio"
,"sigma2_alpha","sigma2_res"
)


library(corrplot)  # For correlation plot

# Convert MCMC list to matrix
samples_mat <- as.matrix(samples) 

# Get parameter names from samples
param_names <- colnames(samples_mat)

# Find parameters that match any in `parToPlot` (including multi-dimensional ones)
matching_params <- param_names[grepl(paste0("^(", paste(parToPlot, collapse = "|"), ")\\[?\\d*\\]?$"), param_names)]

# Check if we found matching parameters
if (length(matching_params) == 0) {
  stop("No matching parameters found. Check parameter names!")
}

# Extract only the relevant parameters
filtered_samples <- samples_mat[, matching_params, drop = FALSE]


# Compute correlation matrix
cor_matrix <- cor(filtered_samples)

# Plot the correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", number.cex = 0.7)
