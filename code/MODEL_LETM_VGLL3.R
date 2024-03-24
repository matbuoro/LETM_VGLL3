LETM<-function(){
  
  for (s in 1:2){
    for (gen in 1:3){
      mu_X[gen,s] ~ dnorm(0,0.001)
      
      sigma_X[gen,s]~dgamma(2,1/s)#dunif(0,10)
      sigma2_X[gen,s] <- pow(sigma_X[gen,s],2)
    }}
  
  s~dchisqr(2)
  
  #======== LIKELIHOOD ========#
  for (i in 1:N) {
    
    ## Environmental cue (X):
    X[i]~dnorm(mu_X[g[i],sex[i]],  1/sigma2_X[g[i],sex[i]]) # Growth distributions
    #X[i]~dnorm(mu_X[g[i], sex[i]], tau_X[g[i], sex[i]]) # Growth distributions
    #X[i]~dnorm(mu[i,1],  1/sigma2_X) # Growth distributions
    #log(mu[i,1]) <- beta[g[i], sex[i], period[i]]
    #X.scaled[i] <- (X[i]-mean(X[1:N]))/sd(X[1:N]) # scaled growth
    #X.scaled[i] <- (X[i]-mean.X[g[i],sex[i]])/sd.X[g[i],sex[i]] # scaled growth
    X.scaled[i] <- (X[i]-mu_X[g[i],sex[i]])/sigma_X[g[i],sex[i]] # scaled growth
    
    ## Phenotypes:
    Y[i]~dbern(p_mat[i])
    p_mat[i]<-phi(z[i]) # probit link function
    #z[i]<-((X.scaled[i] - theta[i]) / sqrt(sigma2_eta[sex[i]]))
    # z[i]<-(((X[i]-mu_X) / sqrt(sigma2_eta[sex[i]] + sigma2_X)) - theta[i]) / sqrt(sigma2_eta[sex[i]]/(sigma2_eta[sex[i]] + sigma2_X))
    # z[i]<-(((X[i]-mu_X[g[i],sex[i]]) / sqrt(sigma2_eta[g[i],sex[i]] + sigma2_X[g[i],sex[i]])) - theta[i]) / sqrt(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]] + sigma2_X[g[i],sex[i]]))
    #z[i]<-(((X[i]-mu_X[g[i],sex[i]]) / sqrt(sigma2_eta[g[i],sex[i]] + 1)) - theta[i]) / sqrt(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]] + 1))
    z[i]<-(((X.scaled[i]) / sqrt(sigma2_eta[g[i],sex[i]] + 1)) - theta[i]) / sqrt(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]] + 1))
    
    ## Threshold
    theta[i]~dnorm(mu[i,2], 1/sigma2_res[g[i],sex[i]])
    #mu[i,2] <- alpha[g[i],sex[i]]
    mu[i,2] <- mu_theta[sex[i]] + alpha[g[i],sex[i]]
    #mu[i,2] <- mu_theta[period[i],sex[i]] + alpha[g[i],sex[i]]
    
    ## Proximate cue (eta)
    lower[i] <- ifelse(Y[i]==1,theta[i],-100)
    upper[i] <- ifelse(Y[i]==0,theta[i], 100)
    mu_eta[i] <- X.scaled[i] / sqrt(sigma2_eta[g[i],sex[i]] + 1) # normalized
    eta[i] ~dnorm(mu_eta[i], 1/(sigma2_eta[g[i],sex[i]]/ (sigma2_eta[g[i],sex[i]]+ 1)));T(lower[i],upper[i]) # /!\ using model as R function, the truncated normal is not a valid R expression. You can fool the R interpreter by inserting a semicolon
    #eta[i] ~dnorm(X[i]-mu_X[g[i],sex[i]]/sqrt(sigma2_eta[g[i],sex[i]] + sigma2_X[g[i],sex[i]]), 1/(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]]+ sigma2_X[g[i],sex[i]])));T(lower[i],upper[i]) # /!\ using model as R function, the truncated normal is not a valid R expression. You can fool the R interpreter by inserting a semicolon
    #eta[i]~dnorm(X[i], 1/sigma2_eta)
    
  } # End of loop i
  
  #mu_theta[1]~dnorm(0, 1) # male
  #mu_theta[2]~dnorm(0, 1)#<- mu_theta[1] + delta_theta
  #delta_theta <- mu_theta[1] - mu_theta[2]#~dnorm(0, 0.001) # female
  mu_theta[1,1]~dnorm(0, 1) # male before 2005
  mu_theta[2,1]~dnorm(0, 1)# male after 2005
  mu_theta[1,2]~dnorm(0, 1) # female before 2005
  mu_theta[2,2]~dnorm(0, 1)# female after 2005
  delta_theta[1] <- mu_theta[1,1] - mu_theta[2,1]#~dnorm(0, 0.001) # female
  delta_theta[2] <- mu_theta[1,2] - mu_theta[2,2]#~dnorm(0, 0.001) # female
  
  for (s in 1:2){    # sex
    # Genotype x Sex effects
    alpha[1,s] <- -a[s] # EE
    alpha[2,s] <- a[s]*k[s] # EL
    alpha[3,s] <- a[s] # LL
    # alpha[1,s]~dnorm(0, 1)
    # alpha[2,s]~dnorm(0, 1)
    # alpha[3,s] <- -alpha[1,s]
    
    # dominance deviation
    #k[s] <- alpha[1,s]/alpha[2,s] # EL
    #d[s] <- k[s]*a # dominance 
    d[s] <- k[s]*a[s] # dominance 
    #k[s]~dnorm(0, 1) # scaled dominance 
    
  } # end loop s
  
  # genotypic values
  a[1]~dnorm(0, 1)
  a[2]~dnorm(0, 1)#<- a[1] + delta_a
  #a[2] <- a[1]
  delta_a <- a[1]-a[2]
  #a[1]<- alpha[3,1]
  #a[2]<- alpha[3,2]
  
  # dominance deviation
  k[1]~dnorm(0, 1) # scaled dominance 
  #k.prior~dnorm(0, .1) # scaled dominance 
  k[2]~dnorm(0, 1)#<- k[1] + delta_k
  delta_k <- k[1] - k[2] 
  
  
  # residual variance of the latent threshold
  for (s in 1:2){
    
    #sigma2_e[s] <- sigma2_eta[s]/(sigma2_eta[s] + sigma2_X)
    
    # Sex variance
    p_sex[s] <- sum(n[1:3,s])/N # female
    #sigma2_sex[s] <- (pow(mu_theta[s] - (mu_theta[1]+p_sex[2]*delta_theta),2))*p_sex[s] # variance sex
   
    mu_alpha[s] <- (a[s]*(-n[1,s]+k[s]*n[2,s]+n[3,s]))/sum(n[1:3,s])
    #mu_alpha[2] <- (a[2]*(-n[1,2]+k[1]*n[2,2]+n[3,2]))/sum(n[1:3,2])
    for (j in 1:3){
      
      # Genotype (VGLL3) x Sex variance
      sigma2_alpha[j,s] <- pow(alpha[j,s] - mu_alpha[s],2)*(n[j,s]/sum(n[1:3,s])) # variance genotype
      
      # Residual genetic variance
      #sigma2_res[j,s] <- sigma_res[j,s]*sigma_res[j,s]
      #sigma_res[j,s] <- chSq[s]/(nu+1)
      #sigma2_res[j,s] <- sigma2_alpha[j,s]*(ratio[2]/(1-ratio[2]))
      
      # TOTAL thresholds variance
      #sigma2_THETA[j,s] <- sigma2_alpha[j,s] + sigma2_res[j,s]#*(n[j,s]/sum(n[1:3,s])) # variance genetique ponderee = variance genotype/sex + residual variance
      
    } # end loop j
    #chSq[s]~dchisqr(nu+1) # male
    # sigma2_ALPHA[s] <- sum(sigma2_alpha[1:3,s])
    # sigma2_RES[s] <- sigma2_res[1,s]*(n[1,s]/sum(n[1:3,s])) + sigma2_res[2,s]*(n[2,s]/sum(n[1:3,s])) + sigma2_res[3,s]*(n[3,s]/sum(n[1:3,s]))
  } # end loop s
  #nu <- 2
  #chSq[1]~dchisqr(nu+1) # male
  #chSq[2] <- chSq[1]#*exp(delta_res) #female
  #delta_res ~dnorm(0,0.001);T(-3,3)
  #chSq.prior~dchisqr(nu+1)
  
  
  for (j in 1:3){
    
    # Variance proximate cue
    #sigma2_eta[j,1] <- sigma2_X[j,1] *(ratio[1]/(1-ratio[1]))
    sigma2_eta[j,1] <- 1*(ratio[1]/(1-ratio[1])) # sigma2_Xscaled = 1
    sigma2_eta[j,2] <- 1*(ratio[2]/(1-ratio[2])) # sigma2_Xscaled = 1
    #sigma2_eta[j,2] <- sigma2_eta[j,1]
    
    
    #sigma2_res[j,1] <- sigma2_alpha[j,1]*(ratio[2]/(1-ratio[2]))
    #sigma2_res[j,2] <- sigma2_res[j,1] # same for sex
    
    # Hypothesis: same  residual  varoance by genotype and sex
    sigma2_res[j,1] <- sum(sigma2_alpha[,1])*(ratio[3]/(1-ratio[3]))
    sigma2_res[j,2] <- sum(sigma2_alpha[,2])*(ratio[4]/(1-ratio[4]))
    #sigma2_res[j,2] <- sigma2_res[j,1] # same for sex
    
  }
  
  
  ratio[1]~dbeta(3,3)
  ratio[2]~dbeta(3,3)
  ratio[3]~dbeta(3,3)
  ratio[4] <- ratio[3]
  #ratio[4]~dbeta(3,3)
  
  #sigma2_eta[1] <- sigma2_X *(ratio[1]/(1-ratio[1]))
  
  #sigma2_eta[1] <- (sigma2_T[1] * (1-h2[1])) - 1 # var(X)=1
  #sigma_eta <- sqrt(sigma2_eta)
  varE[1] <- 1# sigma2_eta[1] + 1#sum(sigma_X[1:3, 1])*(n[j,1]/sum(n[1:3,1])) # Variance environnemental
  #h2[1] ~dbeta(1,1)
  #sigma2_T[1] <- sigma2_G[1] / h2[1]
  sigma2_T[1] <- sigma2_G[1] + varE[1]
  sigma2_G[1] <-  sum(sigma2_alpha[,1]) + sigma2_res[1,1] #sum(sigma2_THETA[1:3,1]) # variance genetic
  h2[1] <- sigma2_G[1] / sigma2_T[1]
  
  # female
  varE[2] <- 1#sigma2_eta[2] + 1
  sigma2_T[2] <- sigma2_G[2] + varE[2]
  sigma2_G[2] <- sum(sigma2_alpha[,2]) + sigma2_res[1,2] #sum(sigma2_THETA[1:3,2]) # variance genetic
  h2[2] <- sigma2_G[2] / sigma2_T[2]
  
  
  
  
  for (s in 1:2){    # sex
    # for (t in 1:2){
    #   p[s,t] <- sum((2*n[1,s,t])+(1*n[2,s,t])+(0*n[3,s,t]))/(2*sum(n[1:3,s,t])) # population frequence allelic for E
    #   q[s,t] <- 1 - p[s,t]
    # }
    p[s] <- sum((2*n[1,s])+(1*n[2,s])+(0*n[3,s]))/(2*sum(n[1:3,s])) # population frequence allelic for E
    q[s] <- 1 - p[s]

    #Average effect of the gene substitution
    gamma[s] <- a[s] + d[s]*(q[s]-p[s])
    #alpha1= q[s]*gamma[s]
    #alpha2= -p[s]*gamma[s]

    #Population mean
    M[s] <- a[s]*(p[s]-q[s])+ (2*d[s]*p[s]*q[s])

    # Phenotypic variance
    # varP = varG + varE = varA + varD + varI + varE
    varT[s] <- varA[s] + varD[s] + varE[s]
    #genetic mean of the population
    varG[s] <- varA[s] + varD[s]
    # additive genetic variance
    varA[s] <- 2*p[s]*q[s]*pow(a[s]+d[s]*(q[s]-p[s]),2)
    #varA[s] <- 2*p*q*pow(a[s]+d[s]*(q-p),2)
    # dominance genetic variance
    varD[s] <- pow(2*p[s]*q[s]*d[s],2)
    #varD[s] <- pow(2*p*q*d[s],2)

    #  h2[s+2] <- varG[s] / varT[s] # total genetic variance contribution = broad-sense heritability
    #  h2[s+4] <- varA[s] / varT[s] # addtive genetic variance contribution = narrow-sense heritability
    #  h2[s+6] <- varD[s] / varT[s] # non-additive genetic variance contribution (dominance + allelic interactions)

  }
  
  # sigma2_TOT <-  sigma2_GENOTYPE + sigma2_ENV # variance total
  # sigma2_ENV <- (varE[1]*p_sex[1]) + (varE[2]*p_sex[2])
  # sigma2_GENOTYPE <- sigma2_GENE + sigma2_SEX
  # #sigma2_ALPHA  <- sigma2_GENE + sigma2_SEX # variance genotype
  # sigma2_SEX <- sigma2_sex[1] + sigma2_sex[2]# variance sex
  # sigma2_VGLL3 <- sum(sigma2_alpha[1:3,1])*p_sex[1] + sum(sigma2_alpha[1:3,2])*p_sex[2]#*n[1:3,1:2])/N) # variance genotype*sex
  # sigma2_RES <- sum(sigma2_res[1:3,1]) * p_sex[1] + sum(sigma2_res[1:3,2])*p_sex[2]
  # #sigma2_RES <- pow(chSq[1]/(nu+1),2) * p_sex[1] + pow(chSq[2]/(nu+1),2)*p_sex[2]
  # sigma2_GENE <- sum(sigma2_THETA[1:3,1]) * p_sex[1] + sum(sigma2_THETA[1:3,2])*p_sex[2]
  # 
  # h[1] <- sigma2_GENOTYPE / sigma2_TOT # total genetic variance contribution = broad-sense heritability
  # h[2] <- (sigma2_ENV - sigma2_eta[1])/ sigma2_TOT
  # h[3] <- sigma2_SEX / sigma2_TOT
  # h[4] <- sigma2_VGLL3 / sigma2_TOT # addtive genetic variance contribution = narrow-sense heritability
  # h[5] <- sigma2_RES / sigma2_TOT # addtive genetic variance contribution = narrow-sense heritability
  # h[6] <- sigma2_GENE / sigma2_TOT # addtive genetic variance contribution = narrow-sense heritability
  # 
  # delta_sigma[1] <- log(sigma2_SEX / sigma2_GENE)
  # delta_sigma[2] <- log(sigma2_SEX / sigma2_VGLL3)
  # delta_sigma[3] <- log(sigma2_SEX / sigma2_RES)
} # end model
