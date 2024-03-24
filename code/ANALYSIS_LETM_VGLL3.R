rm(list=ls())   # Clear memory

# DIRECTORY ####
#setwd("~/Documents/RESEARCH/PROJECTS/VGLL3")

# PACKAGES ####
# Install from source for MacOS, Linux, or Windows:
# install.packages("nimble", repos = "http://r-nimble.org", type = "source")
# the 'type = "source"' is unnecessary for Linx
require(R2jags)
#library(igraph)
library(MCMCvis)
library(writexl)

# DATA ####
source("code/DATA_VGLL3_Scorff.R")
attach(dataToJags)

dataToJags$X.scaled <- NULL


# # Use aggregate and table functions
# result_table=matrix(0,nrow=4,ncol=length(1986:2017)+1); colnames(result_table)<-c("Genotype", paste0(1986:2017)); rownames(result_table)<-c(rep("Male",2),rep("Female",2))
# result_table[,1]<- c("1SW","MSW","1SW","MSW")
# 
# DF <- data.frame(year=dataToJags$year, g=dataToJags$g, sex=dataToJags$sex, Y=dataToJags$Y)
# df <- subset(DF, sex==1 & Y==1)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[1,2:33]<- as.numeric(tmp$count)
# df <- subset(DF, sex==1 & Y==0)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[2,2:33]<- as.numeric(tmp$count)
# 
# df <- subset(DF, sex==2 & Y==1)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[3,2:33]<- as.numeric(tmp$count)
# df <- subset(DF, sex==2 & Y==0)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[4,2:33]<- as.numeric(tmp$count)
# writexl::write_xlsx(as.data.frame(result_table), "results/CountBySexPhenotype.xlsx")
# 
# 
# 
# 
# result_table=matrix(0,nrow=6,ncol=length(1986:2017)+1); colnames(result_table)<-c("Genotype", paste0(1986:2017)); rownames(result_table)<-c(rep("Male",3),rep("Female",3))
# result_table[,1]<- c("EE","EL","LL","EE","EL","LL")
# 
# DF <- data.frame(year=dataToJags$year, g=dataToJags$g, sex=dataToJags$sex, Y=dataToJags$Y)
# df <- subset(DF, sex==1 & g==1)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[1,2:33]<- as.numeric(tmp$count)
# df <- subset(DF, sex==1 & g==2)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[2,2:33]<- as.numeric(tmp$count)
# df <- subset(DF, sex==1 & g==3)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[3,2:33]<- as.numeric(tmp$count)
# df <- subset(DF, sex==2 & g==1)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[4,2:33]<- a &s.numeric(tmp$count)
# df <- subset(DF, sex==2 & g==2)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[5,2:33]<- as.numeric(tmp$count)
# df <- subset(DF, sex==2 & g==3)
# tmp <- aggregate(rep(1,nrow(df)) ~ year, data = df, FUN = sum)
# tmp<-merge(data.frame(year=1986:2017),tmp, by = "year", all = TRUE);colnames(tmp) <- c("year","count")#tmp$`rep(1, nrow(df))`
# result_table[6,2:33]<- as.numeric(tmp$count)
# writexl::write_xlsx(as.data.frame(result_table), "results/CountBySexGenotype.xlsx")


data_empty <- list(N=dataToJags$N, X= dataToJags$X, n=dataToJags$n, g=dataToJags$g, sex=dataToJags$sex)

# MODEL ####
source("code/MODEL_LETM_VGLL3.R")


# ANALYSIS ####
## Parameters ####
parameters <- c(
  #"Y"
  "mu_X","sigma2_X","s"
  #,"beta","eps_X"
  ,"mu_theta", "mu_alpha"#,"eps_theta"
  ,"delta_theta","delta_X","delta_a","delta_k","delta_res"
  ,"alpha", "d","a","k","k.prior"
  ,"chSq","chSq.prior"
  ,"sigma2_eta","sigma2_e","ratio"
  #, "sigma2_theta"
  ,"sigma2_alpha","sigma2_res","sigma2_THETA","h2","sigma2_T","sigma2_G"
  ,"theta","eta"
  #,"gamma","M"
  #,"gap"
  #,"delta","delta_theta","delta_alpha"
  #,"mu_alpha"
  #,"var_alpha","var_theta"
  #,"var_alpha_sex","var_alpha_gen"
  #,"sigma2_GENOTYPE", "sigma2_VGLL3","sigma2_SEX","sigma2_GENE","sigma2_TOT","sigma2_RES","h"
  ,"varG","varA","varD","varT","varE"
  ,"delta_sigma"
  # ,"p"
  
  #,"h"
  #,"psi","lpsi","psi.prior","lpsi.prior"
  #,"theta.pred"
  
  #,"delta_X_sex","delta_sdX_sex","delta_X","delta_sdX", "P_delta","smd_X"
)


## Initial values ####
inits <- function(){
  list(
    #eta=eta, theta=theta
    #MU_theta=1, sigma_theta=2
    # chSq=c(1,NA),delta_res=0
    #, h2=c(0.3,NA, rep(NA,6))
    #mu_theta=c(0,1)#,delta_theta=3.8
    mu_theta=matrix(c(0,0,1,1),nrow=2,ncol=2)
    # ,ratio=c(0.5,0.1)
    #delta_theta=delta_theta,
    , a=c(3,1)#,delta_a=0
    , k=c(-1,0)#,delta_k=1
    #, mu_X=2.7
    #,beta=array(0,dim=c(3,2,2)),sigma_beta=0.1
    #,sigma_X=matrix(.2,3,2)
    ,ratio=rep(0.5, 4)
  )
}

# JAGS ####
samples <- jags.parallel(data=dataToJags,  
                         model.file = LETM,
                         parameters.to.save = parameters,  
                         n.chains = 2,  # Number of chains to run.
                         inits = inits,  # initial values for hyperparameters
                         n.iter = 10000*10,    # n.store * n.thin
                         n.burnin = 1000,   # discard first X iterations
                         n.thin = 10,
                         n.cluster= 2
) # seep every X iterations

save(samples, file="results/RESULTS_vgll3_scorff.RData")

write.csv2(samples$BUGSoutput$summary, file="results/Summary.csv")


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
          filename="results/MCMC_jags_all.pdf",
          wd = "results/",
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE, # add Rhat
          n.eff = TRUE, # add eff sample size
          params = parameters,
          excl=c("theta","eta","Y")
)



# summary(samples)
# gelman.diag(samples)
#
# # MCMCsamples <- as.matrix(samples)
# # plot(samples[ , 'h2'], type = 'l', xlab = 'iteration',  ylab = expression(h^2))
# # plot(density(MCMCsamples[,'h2']))
#
# # # or to use some plots in coda
# # ## Pairs.panel
# panel.hist <- function(x, ...)
# {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(usr[1:2], 0, 1.5) )
#   h <- hist(x, plot = FALSE)
#   breaks <- h$breaks; nB <- length(breaks)
#   y <- h$counts; y <- y/max(y)
#   rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
# }

# pairs(as.matrix(samples),lower.panel=panel.smooth,diag.panel=panel.hist,pch='.')

