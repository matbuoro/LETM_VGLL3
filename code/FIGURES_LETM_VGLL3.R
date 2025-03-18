rm(list=ls())   # Clear memory

# LIBRARY ####
library(mcmcplots)


#____________________DATA______________________#
#source("code/DATA_VGLL3_Scorff.R")
#attach(data)


#____________________CONFIG______________________#
require("NatParksPalettes")
colors <- natparks.pals(name="Banff",n=3,type="discrete")
#colors
# c1 <- rgb(173,216,230, max = 255, alpha = 80, names = "lt.blue")
# c2 <- rgb(211,211,211, max = 255, alpha = 80, names = "lt.grey")
# colors <- c("#00AFBB", "#E7B800", "#FC4E07")
#cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}


# Adds a highest density ellipse to an existing plot
# xy: A matrix or data frame with two columns.
#     If you have to variables just cbind(x, y) them.
# coverage: The percentage of points the ellipse should cover
# border: The color of the border of the ellipse, NA = no border
# fill: The filling color of the ellipse, NA = no fill
# ... : Passed on to the polygon() function
add_hd_ellipse <- function(xy, coverage, border = "blue", fill = NA, ...) {
  library(MASS)
  library(cluster)
  fit <- cov.mve(xy, quantile.used = round(nrow(xy) * coverage))
  points_in_ellipse <- xy[fit$best, ]
  ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
  polygon(ellipse_boundary, border=border, col = fill, ...)
}


library(showtext)
showtext_auto()

female = intToUtf8(9792)
male = intToUtf8(9794)


gene.name=c("EE","EL","LL")
col.sex<-c("lightgrey","darkgrey")



nimble=TRUE
jags=FALSE


# ---- DISTRIBUTIONS -----
if(nimble){
  
  source("code/DATA_VGLL3_Scorff.R")
  attach(data)
  
  X <- data$X
  mu_X <- data$mu_X
  sd_X <- data$sd_X
  nyears <- length(unique(data$year))
  
  
  
  load("results/RESULTS_vgll3_scorff_nimble.RData")
  #stats <- MCMCpstr(samples, func = function(x) quantile(x, probs = c(0.025,0.5 ,0.975)))
  mcmc <- do.call(rbind,samples)
  
  ## pool theta samples between the three chains
  mu_theta <- mcmc[,grep("mu_theta", colnames(mcmc))]
  mu_alphaMale <- mcmc[,"mu_theta[1]"]#samples$BUGSoutput$sims.list$mu_theta[,1]
  mu_alphaFemale <- mcmc[,"mu_theta[2]"]#samples$BUGSoutput$sims.list$mu_theta[,2]
  #mu_alphaMale <- samples$BUGSoutput$sims.list$mu_theta[,1,1]
  #mu_alphaFemale <- samples$BUGSoutput$sims.list$mu_theta[,1,2]
  alpha1m <- mcmc[,"alpha[1, 1]"]
  alpha2m <- mcmc[,"alpha[2, 1]"]
  alpha3m <- mcmc[,"alpha[3, 1]"]
  
  alpha1f <- mcmc[,"alpha[1, 2]"]
  alpha2f <- mcmc[,"alpha[2, 2]"]
  alpha3f <- mcmc[,"alpha[3, 2]"]
  
  theta <- array(,dim=c(length(mu_alphaFemale),3,2))
  theta[,1,1] <- mu_alphaMale + alpha1m
  theta[,2,1] <- mu_alphaMale + alpha2m
  theta[,3,1] <- mu_alphaMale + alpha3m
  
  theta[,1,2] <- mu_alphaFemale + alpha1f
  theta[,2,2] <- mu_alphaFemale + alpha2f
  theta[,3,2] <- mu_alphaFemale + alpha3f
  
  eta <- mcmc[,c(paste0("eta[",1:data$N,"]"))]
  etas <- mcmc[,c(paste0("eta[",1:data$N,"]"))]
  thetas <- mcmc[,c(paste0("theta[",1:data$N,"]"))]
  
  thetas_means <- apply(thetas,2,mean)
  etas_means <- apply(etas,2,mean)
  
  sigma2_alpha <- mcmc[,grep("sigma2_alpha", colnames(mcmc))]
  sigma2_eta <- mcmc[,grep("sigma2_eta", colnames(mcmc))]
  sigma2_res <- mcmc[,grep("sigma2_res", colnames(mcmc))]
  
  # Creating a 3D array sigma2_res
  sigma2_res <- array(sigma2_res, dim = c(nrow(sigma2_res), 3, 2))
  sigma2_alpha <- array(sigma2_alpha, dim = c(nrow(sigma2_res), 3, 2))
  sigma2_eta <- array(sigma2_eta, dim = c(nrow(sigma2_res), 3, 2))
  
  # mu_X <- mcmc[,grep("mu_X", colnames(mcmc))]
  # quantiles_mu_X <- as.data.frame(t(apply(mu_X, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  # colnames(quantiles_mu_X) <- c("2.5%", "25%", "50%", "75%", "97.5%")  
  #  
  # sd_X <- mcmc[,grep("sd_X", colnames(mcmc))]
  # quantiles_sd_X <- as.data.frame(t(apply(mu_X, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  # colnames(quantiles_sd_X) <- c("2.5%", "25%", "50%", "75%", "97.5%")  
  
  ratio <- mcmc[,grep("ratio", colnames(mcmc))]
  a <- mcmc[,c("a[1]","a[2]")]
  k <- mcmc[,c("k[1]","k[2]")]
  h2 <- mcmc[,c("h2[1]","h2[2]")]
  
}

if(jags){
  load("results/RESULTS_vgll3_scorff_jags.RData")
  mcmc <- samples$BUGSoutput$sims.matrix
  mu_alphaMale <- mcmc[,"mu_theta[1]"]#samples$BUGSoutput$sims.list$mu_theta[,1]
  mu_alphaFemale <- mcmc[,"mu_theta[2]"]#samples$BUGSoutput$sims.list$mu_theta[,2]
  #mu_alphaMale <- samples$BUGSoutput$sims.list$mu_theta[,1,1]
  #mu_alphaFemale <- samples$BUGSoutput$sims.list$mu_theta[,1,2]
  alpha1m <- mcmc[,"alpha[1,1]"]
  alpha2m <- mcmc[,"alpha[2,1]"]
  alpha3m <- mcmc[,"alpha[3,1]"]
  
  alpha1f <- mcmc[,"alpha[1,2]"]
  alpha2f <- mcmc[,"alpha[2,2]"]
  alpha3f <- mcmc[,"alpha[3,2]"]

  theta <- array(,dim=c(length(mu_alphaFemale),3,2))
  theta[,1,1] <- mu_alphaMale + alpha1m
  theta[,2,1] <- mu_alphaMale + alpha2m
  theta[,3,1] <- mu_alphaMale + alpha3m
  
  theta[,1,2] <- mu_alphaFemale + alpha1f
  theta[,2,2] <- mu_alphaFemale + alpha2f
  theta[,3,2] <- mu_alphaFemale + alpha3f
  
  eta <- samples$BUGSoutput$sims.list$eta
  thetas <- samples$BUGSoutput$sims.list$theta
  etas <- (samples$BUGSoutput$sims.list$eta)
  
  thetas_means <- samples$BUGSoutput$mean$theta
  etas_means <- samples$BUGSoutput$mean$eta
  
  # alphas <- stats$alpha
  # mu_theta <- samples$BUGSoutput$median$mu_theta
  # 
  # sigma2_eta <- stats$sigma2_eta
  # sigma2_theta <- stats$sigma2_theta
  # sigma2_THETA <- stats$sigma2_THETA
  # 
  sigma2_THETA <- samples$BUGSoutput$sims.list$sigma2_THETA
  sigma2_alpha <- samples$BUGSoutput$sims.list$sigma2_alpha
  sigma2_eta <- samples$BUGSoutput$sims.list$sigma2_eta
  sigma2_res <- samples$BUGSoutput$sims.list$sigma2_res
  
  ratio <- samples$BUGSoutput$sims.list$ratio 
  a <- samples$BUGSoutput$sims.list$a 
  k <- samples$BUGSoutput$sims.list$k
  h2 <- samples$BUGSoutput$sims.list$h2
}


if(nimble){
pdf("results/FIGURES_VGLL3_Scorff_nimble.pdf")
}else{
pdf("results/FIGURES_VGLL3_Scorff_jags.pdf")
}



par(mfrow=c(2,2))
# male
plot(sigma2_res[,1,1], sigma2_eta[,1,1],pch=".",xlab="Variance Gen residual", ylab="Variance proximate cue", main="Male")
plot(sigma2_alpha[,1,1], sigma2_eta[,1,1],pch=".",xlab="Variance Vgll3", ylab="Variance proximate cue", main="Male")
#plot(sigma2_THETA[,1,1], sigma2_eta[,1,1],pch=".",xlab="Variance Gen tot", ylab="Variance proximate cue")

# female
plot(sigma2_res[,1,2], sigma2_eta[,1,2],pch=".",xlab="Variance Gen residual", ylab="Variance proximate cue", main="Female")
plot(sigma2_alpha[,1,2], sigma2_eta[,1,2],pch=".",xlab="Variance Vgll3", ylab="Variance proximate cue", main="Female")
#plot(sigma2_THETA[,1,2], sigma2_eta[,1,2],pch=".",xlab="Variance Gen tot", ylab="Variance proximate cue")


par(mfrow=c(1,3))
plot(sigma2_res[,1,1], sigma2_res[,1,2],pch=".",xlab="Variance Gen residual male", ylab="ariance Gen residual female")
plot(sigma2_res[,1,1], sigma2_alpha[,1,1],pch=".",xlab="Variance Gen residual", ylab="Variance VGLL3")
plot(sigma2_res[,1,1], sigma2_eta[,1,1],pch=".",xlab="Variance Gen residual", ylab="Variance eta")


par(mfcol=c(1,3))
#par(mfcol=c(1,3),
#oma=c(1,1,0,0),
#mar=c(1,1,1,0),
#mgp=c(3,1,0)
#xpd=NA
#)
ylim=range(new.df$X)
xlim=c(-1.5,1.5)

for (gene in 1:3){
  
  side=c(-1,1)
  #Make the plot
  if(gene==1){
    plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="Marine growth", xaxt="n")#,main=gene.name[gene])
    legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex, title ="Growth distributions")
    axis(1, line=0, at=0,labels=c("EE"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
  } 
  if(gene==2){
    plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n",yaxt="n")#,main=gene.name[gene])

    axis(1, line=0, at=0,labels=c("EL"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
    axis(1, line=2, at=0,labels=c("GENOTYPE"),tick = FALSE,cex.axis=2)#, "EL", "LL"))
  } 
  
  if(gene==3){
    plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n", yaxt="n")#,main=gene.name[gene])
    axis(1, line=0, at=0,labels=c("LL"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
  }
  #axis(side = 1, at=0,labels = type)
  #axis(2, at=1:length(seq(40,140,10)),labels=seq(40,140,10), col.axis="red", las=2)
  
  text(-0.5,1,male,cex=2)
  text(0.5,1,female,cex=2)
  
  
  for (sexe in 1:2){
    
    tmp <- X[g==gene & sex==sexe]

    polygon(density(tmp)$y * side[sexe],density(tmp)$x , main="" , xlab="",ylim=ylim , xaxt="n", col=t_col(col.sex[sexe],20) , border=NA)
    segments(0,median(tmp),max(density(tmp)$y) * side[sexe],median(tmp),col="#99999950")
    points(0,mean(tmp),pch=16,col=col.sex[sexe],cex=3)
    points(rnorm(length(tmp),0.5*side[sexe],0.025), tmp,pch=16,col=t_col(col.sex[sexe],10) )
    
  }
}


#id <-which(data$X>550)



### ETA
### ## TODO




par(mfrow=c(2,2))
#boxplot(X~g + sex, data=data)
plot(NULL, xlim=c(0,4),ylim=range(data$X),xlab="Genotype", ylab="Size",xaxt="n")
axis(1, at=1:3,labels=c("EE", "EL", "LL"))
legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")

for (j in 1:3){
  # MALE
  segments(j-0.1, quantile(data$X[g==j & sex==1], probs = 0.025),j-0.1, quantile(data$X[g==j & sex==1], probs = 0.975),col=colors[j])
  segments(j-0.1, quantile(data$X[g==j & sex==1], probs = 0.25),j-0.1, quantile(data$X[g==j & sex==1], probs = 0.75),col=colors[j],lwd=2)
  points(j-0.1,median(data$X[g==j & sex==1]), pch=24, col=colors[j],bg="white")
  text(j-0.1,quantile(data$X[g==j & sex==1], probs = 0.995),male,col=colors[j])

  # FEMALE
  segments(j+0.1, quantile(data$X[g==j & sex==2], probs = 0.025),j+0.1, quantile(data$X[g==j & sex==2], probs = 0.975),col=colors[j])
  segments(j+0.1, quantile(data$X[g==j & sex==2], probs = 0.25),j+0.1, quantile(data$X[g==j & sex==2], probs = 0.75),col=colors[j],lwd=2)
  points(j+0.1,median(data$X[g==j & sex==2]), pch=21, col=colors[j],bg="white")
  text(j+0.1,quantile(data$X[g==j & sex==2], probs = 0.995),female,col=colors[j])
}



#par(mfrow=c(2,2))
plot(NULL, xlim=c(0,4),ylim=c(-4.5,4.5),xlab="Genotype", ylab="Means of maturation thresholds",xaxt="n")
axis(1, at=1:3,labels=c("EE", "EL", "LL"))
legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
#abline(h=mu_alphas[1,"50%"],lty=3,col="darkgrey");text(4,mu_alphas[1,"50%"]+1,male,col="darkgrey")
#abline(h=mu_alphas[2,"50%"],lty=5,col="darkgrey");text(4,mu_alphas[2,"50%"]+1,female,col="darkgrey")
for (j in 1:3){
  #MALE
  segments(1-0.1, apply(theta[,,1], 2, quantile, probs=0.5)[1],3-0.1, apply(theta[,,1], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
  segments(j-0.1, apply(theta[,,1], 2, quantile, probs=0.025)[j],j-0.1, apply(theta[,,1], 2, quantile, probs=0.975)[j],col=colors[j])
  segments(j-0.1, apply(theta[,,1], 2, quantile, probs=0.25)[j],j-0.1, apply(theta[,,1], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
  points(j-0.1,apply(theta[,,1], 2, quantile, probs=0.5)[j], pch=24, col=colors[j],bg="white")
  #text(j-0.1,apply(theta[,,1], 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
  
  
  #FEMALE
  segments(1+0.1, apply(theta[,,2], 2, quantile, probs=0.5)[1],3+0.1, apply(theta[,,2], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
  segments(j+0.1, apply(theta[,,2], 2, quantile, probs=0.025)[j],j+0.1, apply(theta[,,2], 2, quantile, probs=0.975)[j],col=colors[j])
  segments(j+0.1, apply(theta[,,2], 2, quantile, probs=0.25)[j],j+0.1, apply(theta[,,2], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
  points(j+0.1,apply(theta[,,2], 2, quantile, probs=0.5)[j], pch=21, col=colors[j],bg="white")
  #text(j+0.1,apply(theta[,,2], 2, quantile, probs=0.99)[j]+1,female,col=colors[j])
  
}




## REACTION NORMS
compute_reaction_norms <- function(mcmc, X, sigma2_eta, theta) {
  # Define X.sim and X.obs
  X.sim <- seq(200, 600, 10)
  X.pred <- scale(X.sim, center = mean(X), scale = sd(X))
  X.obs <- seq(min(X), max(X), 10)
  X.obsc <- scale(X.obs, center = mean(X), scale = sd(X))
  
  # Initialize arrays
  p.pred <- array(NA, dim = c(length(X.pred), 3, 2))
  p.obs <- array(NA, dim = c(length(X.obs), 3, 2))
  
  if(nimble){
    # Compute quantiles for mu_X
    # mu_X <- mcmc[, grep("mu_X", colnames(mcmc))]
    # quantiles_mu_X <- apply(mu_X, 2, quantile, probs = c(0.05, 0.5, 0.95))
    # mu_X <- matrix(quantiles_mu_X[2,], nrow = 3, ncol = 2, byrow = FALSE)
    # 
    # # Compute quantiles for sd_X
    # sd_X <- sqrt(mcmc[, grep("sigma2_X", colnames(mcmc))])
    # quantiles_sd_X <- apply(sd_X, 2, quantile, probs = c(0.05, 0.5, 0.95))
    # sd_X <- matrix(quantiles_sd_X[2,], nrow = 3, ncol = 2, byrow = FALSE)
    
    # Compute median sigma2_eta
    sigma2_eta <- mcmc[,grep("sigma2_eta", colnames(mcmc))]
    sigma2_eta_med <- apply(sigma2_eta, 2, quantile, probs = c(0.05, 0.5, 0.95))
  }else{
    #mu_X <- samples$BUGSoutput$median$mu_X
    #sd_X <- sqrt(samples$BUGSoutput$median$sigma2_X)
    sigma2_eta <- samples$BUGSoutput$median$sigma2_eta
  }
  
  for (i in 1:3) {
    for (j in 1:2) {
      
      theta_med <- apply(theta[,,j], 2, quantile, probs = c(0.05, 0.5, 0.95))
      sigma2_eta_idx <- i + (j - 1) * 3
      
      X.scaled <- (X.sim - mu_X[i, j]) / sd_X[i, j]
      z_pred <- ((X.scaled / sqrt(sigma2_eta_med[2, sigma2_eta_idx] + 1)) - theta_med[2, i]) /
        sqrt(sigma2_eta_med[2, sigma2_eta_idx] / (sigma2_eta_med[2, sigma2_eta_idx] + 1))
      
      X.scaled <- (X.obs - mu_X[i, j]) / sd_X[i, j]
      z_obs <- ((X.scaled / sqrt(sigma2_eta_med[2, sigma2_eta_idx] + 1)) - theta_med[2, i]) /
        sqrt(sigma2_eta_med[2, sigma2_eta_idx] / (sigma2_eta_med[2, sigma2_eta_idx] + 1))
      
      #if (j == 1) {
      p.pred[, i, j] <- pnorm(z_pred)
      #} else {
      p.obs[, i, j] <- pnorm(z_obs)
    }
  }
  
  
  return(list(p.pred = p.pred, p.obs = p.obs, X.sim=X.sim, X.pred=X.pred, X.obsc=X.obsc))
}

result <- compute_reaction_norms(mcmc, X, sigma2_eta, theta)

X.sim <- result$X.sim
p.pred <- result$p.pred
X.pred <- result$X.pred
X.obsc <- result$X.obsc
p.obs <- result$p.obs


#par(mfrow=c(1,2))

hist((data$X[data$sex==1] - mean(data$X))/sd(data$X)
     , xlim=range(X.pred)
     ,border=t_col("lightgrey",20), lty="blank"
     ,xlab="",ylab="", xaxt='n', yaxt='n',main=""
     , freq = FALSE,col = t_col("lightgrey",20)) # Save 2nd histogram data
#plot(hgB, col = c2)

# prepare graphics to add second plot
par(new = TRUE)

ylab="Maturation probability at 1SW"
plot(NULL, xlim=range(X.pred),ylim=c(0,1)
     , xaxt='n'
     ,ylab=ylab, xlab="Growth", main=male)
axis(1,at=X.pred,labels=X.sim)
legend("topleft", legend=c("EE", "EL", "LL"),fill=colors, bty="n",border = colors)
lines(X.pred, p.pred[,1,1],lty=1,lwd=1,col=colors[1])
lines(X.obsc, p.obs[,1,1],lty=1,lwd=3,col=colors[1])

lines(X.pred, p.pred[,2,1],lty=1,lwd=1,col=colors[2])
lines(X.obsc, p.obs[,2,1],lty=1,lwd=3,col=colors[2])

lines(X.pred, p.pred[,3,1],lty=1,lwd=1,col=colors[3])
lines(X.obsc, p.obs[,3,1],lty=1,lwd=3,col=colors[3])



hist((data$X[sex==2] - mean(data$X))/sd(data$X)
     , xlim=range(X.pred)
     ,border=t_col("lightgrey",20),lty="blank"
     ,xlab="",ylab="", xaxt='n', yaxt='n',main=""
     , freq = FALSE,col = t_col("lightgrey",20)) # Save 2nd histogram data
#plot(hgB, col = c2, add = TRUE)
# prepare graphics to add second plot
par(new = TRUE)
plot(NULL, xlim=range(X.pred),ylim=c(0,1),ylab=ylab
     , xaxt='n'
     , xlab="Growth", main=female)
axis(1,at=X.pred,labels=X.sim)
legend("topleft", legend=c("EE", "EL", "LL"),fill=colors, bty="n",border = colors)
lines(X.pred, p.pred[,1,2],lty=1,lwd=1,col=colors[1])
lines(X.obsc, p.obs[,1,2],lty=1,lwd=3,col=colors[1])

lines(X.pred, p.pred[,2,2],lty=1,lwd=1,col=colors[2])
lines(X.obsc, p.obs[,2,2],lty=1,lwd=3,col=colors[2])

lines(X.pred, p.pred[,3,2],lty=1,lwd=1,col=colors[3])
lines(X.obsc, p.obs[,3,2],lty=1,lwd=3,col=colors[3])




#### VARIANCES


par(mfrow=c(2,2))

# plot(NULL, xlim=c(0,4),ylim=c(0,2),xlab="Genotype", ylab="Total variances of maturation thresholds",xaxt="n")
# axis(1, at=1:3,labels=c("EE", "EL", "LL"))
# legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
# #abline(h=mu_alphas[1,"50%"],lty=3,col="darkgrey");text(4,mu_alphas[1,"50%"]+1,male,col="darkgrey")
# #abline(h=mu_alphas[2,"50%"],lty=5,col="darkgrey");text(4,mu_alphas[2,"50%"]+1,female,col="darkgrey")
# for (j in 1:3){
#   #MALE
#   #segments(1-0.1, apply(sigma2_THETA[,,1], 2, quantile, probs=0.5)[1],3-0.1, apply(sigma2_THETA[,,1], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
#   segments(j-0.1, apply(sigma2_THETA[,,1], 2, quantile, probs=0.025)[j],j-0.1, apply(sigma2_THETA[,,1], 2, quantile, probs=0.975)[j],col=colors[j])
#   segments(j-0.1, apply(sigma2_THETA[,,1], 2, quantile, probs=0.25)[j],j-0.1, apply(sigma2_THETA[,,1], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
#   points(j-0.1,apply(sigma2_THETA[,,1], 2, quantile, probs=0.5)[j], pch=24, col=colors[j],bg="white")
#   text(j-0.1,apply(sigma2_THETA[,,1], 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
#   
#   
#   #FEMALE
#   #segments(1+0.1, apply(sigma2_THETA[,,2], 2, quantile, probs=0.5)[1],3+0.1, apply(sigma2_THETA[,,2], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
#   segments(j+0.1, apply(sigma2_THETA[,,2], 2, quantile, probs=0.025)[j],j+0.1, apply(sigma2_THETA[,,2], 2, quantile, probs=0.975)[j],col=colors[j])
#   segments(j+0.1, apply(sigma2_THETA[,,2], 2, quantile, probs=0.25)[j],j+0.1, apply(sigma2_THETA[,,2], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
#   points(j+0.1,apply(sigma2_THETA[,,2], 2, quantile, probs=0.5)[j], pch=21, col=colors[j],bg="white")
#   text(j+0.1,apply(sigma2_THETA[,,2], 2, quantile, probs=0.99)[j]+1,female,col=colors[j])
#   
# }




plot(NULL, xlim=c(0,4),ylim=c(0,2),xlab="Genotype", ylab="VGLL3 variances",xaxt="n")
axis(1, at=1:3,labels=c("EE", "EL", "LL"))
legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
#abline(h=mu_alphas[1,"50%"],lty=3,col="darkgrey");text(4,mu_alphas[1,"50%"]+1,male,col="darkgrey")
#abline(h=mu_alphas[2,"50%"],lty=5,col="darkgrey");text(4,mu_alphas[2,"50%"]+1,female,col="darkgrey")
for (j in 1:3){
  #MALE
  #segments(1-0.1, apply(sigma2_alpha[,,1], 2, quantile, probs=0.5)[1],3-0.1, apply(sigma2_alpha[,,1], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
  segments(j-0.1, apply(sigma2_alpha[,,1], 2, quantile, probs=0.025)[j],j-0.1, apply(sigma2_alpha[,,1], 2, quantile, probs=0.975)[j],col=colors[j])
  segments(j-0.1, apply(sigma2_alpha[,,1], 2, quantile, probs=0.25)[j],j-0.1, apply(sigma2_alpha[,,1], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
  points(j-0.1,apply(sigma2_alpha[,,1], 2, quantile, probs=0.5)[j], pch=24, col=colors[j],bg="white")
  text(j-0.1,apply(sigma2_alpha[,,1], 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
  
  
  #FEMALE
  #segments(1+0.1, apply(sigma2_alpha[,,2], 2, quantile, probs=0.5)[1],3+0.1, apply(sigma2_alpha[,,2], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
  segments(j+0.1, apply(sigma2_alpha[,,2], 2, quantile, probs=0.025)[j],j+0.1, apply(sigma2_alpha[,,2], 2, quantile, probs=0.975)[j],col=colors[j])
  segments(j+0.1, apply(sigma2_alpha[,,2], 2, quantile, probs=0.25)[j],j+0.1, apply(sigma2_alpha[,,2], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
  points(j+0.1,apply(sigma2_alpha[,,2], 2, quantile, probs=0.5)[j], pch=21, col=colors[j],bg="white")
  text(j+0.1,apply(sigma2_alpha[,,2], 2, quantile, probs=0.99)[j]+1,female,col=colors[j])
  
}




plot(NULL, xlim=c(0,4),ylim=c(0,20),xlab="Genotype", ylab="Residual genetic variances",xaxt="n")
axis(1, at=1:3,labels=c("EE", "EL", "LL"))
legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
#abline(h=mu_alphas[1,"50%"],lty=3,col="darkgrey");text(4,mu_alphas[1,"50%"]+1,male,col="darkgrey")
#abline(h=mu_alphas[2,"50%"],lty=5,col="darkgrey");text(4,mu_alphas[2,"50%"]+1,female,col="darkgrey")
for (j in 1:3){
  #MALE
  #segments(1-0.1, apply(sigma2_res[,,1], 2, quantile, probs=0.5)[1],3-0.1, apply(sigma2_res[,,1], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
  segments(j-0.1, apply(sigma2_res[,,1], 2, quantile, probs=0.025)[j],j-0.1, apply(sigma2_res[,,1], 2, quantile, probs=0.975)[j],col=colors[j])
  segments(j-0.1, apply(sigma2_res[,,1], 2, quantile, probs=0.25)[j],j-0.1, apply(sigma2_res[,,1], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
  points(j-0.1,apply(sigma2_res[,,1], 2, quantile, probs=0.5)[j], pch=24, col=colors[j],bg="white")
  text(j-0.1,apply(sigma2_res[,,1], 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
  
  
  #FEMALE
  #segments(1+0.1, apply(sigma2_res[,,2], 2, quantile, probs=0.5)[1],3+0.1, apply(sigma2_res[,,2], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
  segments(j+0.1, apply(sigma2_res[,,2], 2, quantile, probs=0.025)[j],j+0.1, apply(sigma2_res[,,2], 2, quantile, probs=0.975)[j],col=colors[j])
  segments(j+0.1, apply(sigma2_res[,,2], 2, quantile, probs=0.25)[j],j+0.1, apply(sigma2_res[,,2], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
  points(j+0.1,apply(sigma2_res[,,2], 2, quantile, probs=0.5)[j], pch=21, col=colors[j],bg="white")
  text(j+0.1,apply(sigma2_res[,,2], 2, quantile, probs=0.99)[j]+1,female,col=colors[j])
  
}



par(mfrow=c(1,3))


if(nimble){ratio <- mcmc[,grep("ratio", colnames(mcmc))]}else{
  ratio <- samples$BUGSoutput$sims.list$ratio 
}


plot(NULL, xlim=c(0,3),ylim=c(0,1),xlab="", ylab="Contribution proximate cue to environmental variance",xaxt="n",main="Environment")
#axis(1, at=c(1.5,3.5),labels=c("Envir", "Genetic"),las=2)
axis(1, at=c(1,2),labels=c("Male", "Female"),las=2)
#axis(1, at=1:4,labels=c("ratio_eta", "ratio_eta", "ratio_res", "ratio_res"),las=2)
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
for (j in 1:2){
  segments(j, apply(ratio, 2, quantile, probs=0.025)[j],j, apply(ratio, 2, quantile, probs=0.975)[j],col="black")
  segments(j, apply(ratio, 2, quantile, probs=0.25)[j],j, apply(ratio, 2, quantile, probs=0.75)[j],col="black",lwd=2)
  points(j,apply(ratio, 2, quantile, probs=0.5)[j], pch=21, col="black",bg="white")
  text(j,apply(ratio, 2, quantile, probs=0.99)[j]+1,male,col="black")
}

plot(NULL, xlim=c(0,3),ylim=c(0,1),xlab="", ylab="Contribution of VGLL3 to genetic variance",xaxt="n",main="Genetic")
axis(1, at=c(1,2),labels=c("Male", "Female"),las=2)
#axis(1, at=1:4,labels=c("ratio_eta", "ratio_eta", "ratio_res", "ratio_res"),las=2)
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
for (j in 3:4){
  segments(j-2, apply(1-ratio, 2, quantile, probs=0.025)[j],j-2, apply(1-ratio, 2, quantile, probs=0.975)[j],col="black")
  segments(j-2, apply(1-ratio, 2, quantile, probs=0.25)[j],j-2, apply(1-ratio, 2, quantile, probs=0.75)[j],col="black",lwd=2)
  points(j-2,apply(1-ratio, 2, quantile, probs=0.5)[j], pch=21, col="black",bg="white")
  text(j-2,apply(1-ratio, 2, quantile, probs=0.99)[j]+1,female,col="black")
}
  
# Contribution total variance = Percentages of the total phenotypic variance
#h2 <- samples$BUGSoutput$sims.list$h2
#h2 <- mcmc[,paste0("h2[",1:2,"]")]
plot(NULL, xlim=c(0,3),ylim=c(0,1),xlab="", ylab="Genetic contribution to total variance (h2)",xaxt="n",main="Heritability")
axis(1, at=1:2
     ,labels=c(male, female))
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
for (j in 1:2){
  segments(j, apply(h2, 2, quantile, probs=0.025)[j],j, apply(h2, 2, quantile, probs=0.975)[j],col=1)
  segments(j, apply(h2, 2, quantile, probs=0.25)[j],j, apply(h2, 2, quantile, probs=0.75)[j],col=1,lwd=2)
  points(j,apply(h2, 2, quantile, probs=0.5)[j], pch=21, col=1,bg="white")
  #text(j-0.1,apply(mu_alphas, 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
}


# #h <- samples$BUGSoutput$sims.list$h
# h <- mcmc[,paste0("h[",c(1,3,4,5),"]")]
# plot(NULL, xlim=c(0,5),ylim=c(0,1),xlab="", ylab="Total variance contribution",xaxt="n")
# axis(1, at=1:4
#      ,labels=c("GENETIC", "SEX","VGLL3","Residual")#,"ENV")
#      ,las=2)
# #legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
# #for (j in 1:5){
# i=0
# for (j in 1:4){
#   i=i+1
#   segments(j, apply(h, 2, quantile, probs=0.025)[j],j, apply(h, 2, quantile, probs=0.975)[j],col=1)
#   segments(j, apply(h, 2, quantile, probs=0.25)[j],j, apply(h, 2, quantile, probs=0.75)[j],col=1,lwd=2)
#   points(j,apply(h, 2, quantile, probs=0.5)[j], pch=21, col=1,bg="white")
#   #text(j-0.1,apply(mu_alphas, 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
# }



# sigma2_eta <- samples$BUGSoutput$sims.list$sigma2_eta
# plot(NULL, xlim=c(0,4),ylim=c(0,15),xlab="Genotype", ylab="Variances of proximate cue",xaxt="n")
# axis(1, at=1:3,labels=c("EE", "EL", "LL"))
# legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
# #abline(h=mu_alphas[1,"50%"],lty=3,col="darkgrey");text(4,mu_alphas[1,"50%"]+1,male,col="darkgrey")
# #abline(h=mu_alphas[2,"50%"],lty=5,col="darkgrey");text(4,mu_alphas[2,"50%"]+1,female,col="darkgrey")
# for (j in 1:3){
#   #MALE
#   #segments(1-0.1, apply(sigma2_eta[,,1], 2, quantile, probs=0.5)[1],3-0.1, apply(sigma2_eta[,,1], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
#   segments(j-0.1, apply(sigma2_eta[,,1], 2, quantile, probs=0.025)[j],j-0.1, apply(sigma2_eta[,,1], 2, quantile, probs=0.975)[j],col=colors[j])
#   segments(j-0.1, apply(sigma2_eta[,,1], 2, quantile, probs=0.25)[j],j-0.1, apply(sigma2_eta[,,1], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
#   points(j-0.1,apply(sigma2_eta[,,1], 2, quantile, probs=0.5)[j], pch=24, col=colors[j],bg="white")
#   text(j-0.1,apply(sigma2_eta[,,1], 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
#   
#   
#   #FEMALE
#   #segments(1+0.1, apply(sigma2_eta[,,2], 2, quantile, probs=0.5)[1],3+0.1, apply(sigma2_eta[,,2], 2, quantile, probs=0.5)[3],col="darkgrey",lty=2)
#   segments(j+0.1, apply(sigma2_eta[,,2], 2, quantile, probs=0.025)[j],j+0.1, apply(sigma2_eta[,,2], 2, quantile, probs=0.975)[j],col=colors[j])
#   segments(j+0.1, apply(sigma2_eta[,,2], 2, quantile, probs=0.25)[j],j+0.1, apply(sigma2_eta[,,2], 2, quantile, probs=0.75)[j],col=colors[j],lwd=2)
#   points(j+0.1,apply(sigma2_eta[,,2], 2, quantile, probs=0.5)[j], pch=21, col=colors[j],bg="white")
#   text(j+0.1,apply(sigma2_eta[,,2], 2, quantile, probs=0.99)[j]+1,female,col=colors[j])
#   
# }



#### GENETIC COMPONENTS
par(mfrow=c(1,3))

mu_alphas <- cbind(mu_alphaMale,mu_alphaFemale)
# Mean threshiolds
#mu_alphas <- samples$BUGSoutput$sims.list$mu_alpha
plot(NULL, xlim=c(0,3),ylim=range(mu_alphas),xlab="Sex", ylab="Mean Maturation thresholds",xaxt="n")
axis(1, at=1:2,labels=c(male, female))
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
#abline(h=mu_alphas[1],lty=3,col="lightgrey");text(4,mu_alphas[1]+1,male,col="lightgrey")
#abline(h=mu_alphas[2],lty=5,col="lightgrey");text(4,mu_alphas[2]+1,female,col="lightgrey")
for (j in 1:2){
  #MALE
  segments(j, apply(mu_alphas, 2, quantile, probs=0.025)[j],j, apply(mu_alphas, 2, quantile, probs=0.975)[j],col=1)
  segments(j, apply(mu_alphas, 2, quantile, probs=0.25)[j],j, apply(mu_alphas, 2, quantile, probs=0.75)[j],col=1,lwd=2)
  points(j,apply(mu_alphas, 2, quantile, probs=0.5)[j], pch=21, col=1,bg="white")
  #text(j,apply(mu_alphas, 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
}


# Additive effect
#am <- c(samples$chain1[,"a[1]"], samples$chain2[,"a[1]"],samples$chain3[,"a[1]"])
#af <- c(samples$chain1[,"a[2]"], samples$chain2[,"a[2]"],samples$chain3[,"a[2]"])
#a <- cbind(am,af)
#a <- samples$BUGSoutput$sims.list$a
plot(NULL, xlim=c(0,3),ylim=c(0,4),xlab="Sex", ylab="Genotypic value",xaxt="n")
axis(1, at=1:2,labels=c(male, female))
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
#abline(h=mu_alphas[1],lty=3,col="lightgrey");text(4,mu_alphas[1]+1,male,col="lightgrey")
#abline(h=mu_alphas[2],lty=5,col="lightgrey");text(4,mu_alphas[2]+1,female,col="lightgrey")
for (j in 1:2){
  #MALE
  segments(j, apply(a, 2, quantile, probs=0.025)[j],j, apply(a, 2, quantile, probs=0.975)[j],col=1)
  segments(j, apply(a, 2, quantile, probs=0.25)[j],j, apply(a, 2, quantile, probs=0.75)[j],col=1,lwd=2)
  points(j,apply(a, 2, quantile, probs=0.5)[j], pch=21, col=1,bg="white")
  #text(j-0.1,apply(mu_alphas, 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
}


# Dominance effect
#k <- samples$BUGSoutput$sims.list$k
# km <- c(samples$chain1[,"k[1]"], samples$chain2[,"k[1]"],samples$chain3[,"k[1]"])
# kf <- c(samples$chain1[,"k[2]"], samples$chain2[,"k[2]"],samples$chain3[,"k[2]"])
# k <- cbind(km,kf)
#k <- samples$BUGSoutput$sims.list$k
plot(NULL, xlim=c(0,3),ylim=range(k),xlab="Sex", ylab="Dominance deviation (scaled)",xaxt="n")
axis(1, at=1:2,labels=c(male, female))
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
abline(h=0,lty=2,col="darkgrey")
#abline(h=mu_alphas[2],lty=5,col="lightgrey");text(4,mu_alphas[2]+1,female,col="lightgrey")
for (j in 1:2){
  #MALE
  segments(j, apply(k, 2, quantile, probs=0.025)[j],j, apply(k, 2, quantile, probs=0.975)[j],col=1)
  segments(j, apply(k, 2, quantile, probs=0.25)[j],j, apply(k, 2, quantile, probs=0.75)[j],col=1,lwd=2)
  points(j,apply(k, 2, quantile, probs=0.5)[j], pch=21, col=1,bg="white")
  #text(j-0.1,apply(mu_alphas, 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
}

#dev.off()














# par(mfrow=c(1,1))
# delta_sigma <- mcmc[,paste0("delta_sigma[",c(1,2),"]")]
# # delta_sigma[1] <- log(sigma2_SEX / sigma2_GENE)
# # delta_sigma[2] <- log(sigma2_SEX / sigma2_VGLL3)
# # delta_sigma[3] <- log(sigma2_SEX / sigma2_RES)
# plot(NULL, xlim=c(0,3),ylim=c(-1,1),xlab="", ylab=" Ratio variance (log)",xaxt="n")
# axis(1, at=1:2
#      ,labels=c("sigma2_SEX / sigma2_GENE", "sigma2_SEX / sigma2_VGLL3")
#      ,las=2)
# #legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
# abline(h=0,lty=5,col="lightgrey")
# for (j in 1:2){
#   segments(j, apply(delta_sigma, 2, quantile, probs=0.025)[j],j, apply(delta_sigma, 2, quantile, probs=0.975)[j],col=1)
#   segments(j, apply(delta_sigma, 2, quantile, probs=0.25)[j],j, apply(delta_sigma, 2, quantile, probs=0.75)[j],col=1,lwd=2)
#   points(j,apply(delta_sigma, 2, quantile, probs=0.5)[j], pch=21, col=1,bg="white")
#   #text(j-0.1,apply(mu_alphas, 2, quantile, probs=0.99)[j]+1,male,col=colors[j])
# }




###### Plot 2
#sigma2_eta <- samples$BUGSoutput$median$sigma2_eta
# sigma2_eta <- apply(mcmc[,paste0("sigma2_eta[",1:2,"]")],2,median)
# sigma_eta <- sqrt(sigma2_eta)
# sex <- constants$sex
# eta=NULL
# for (i in 1:constants$N){
#   eta[i] <- rnorm(1,X.scaled[i], (sigma_eta[sex[i]]))
# }

#eta <- samples$BUGSoutput$median$eta
#eta <- samples$BUGSoutput$sims.list$eta
#X.scaled <- (X-mean(X))/sd(X)



#mu_X_med <- array(quantiles_mu_X[,"50%"],dim=c(3,2))
#sd_X_med <- array(quantiles_sd_X[,"50%"],dim=c(3,2))
X.scaled=NULL
for (i in 1: length(X)){
  X.scaled[i] <- (X[i]-mu_X[data$g[i],data$sex[i]])/sd_X[data$g[i],data$sex[i]]
}


par(mfrow=c(3,2))
for (j in 1:6){
  
# Data frame
smp <- sample(1:nrow(eta),1,replace=TRUE)
df <- data.frame(x = eta[smp,], y = X.scaled, group = sex, genotype = g)

plot(df$x,df$y,pch=df$genotype,col=sex,
     xlim=c(-4,4),ylim=c(-4,4)
     ,xlab="Proximate cue",ylab="Marine growth (scaled)"
     ,main=paste0("Iteration: ", smp)
     )
legend("topleft", legend=c("Male", "Female"),fill=1:2, bty="n",border = NA)
legend("topright", legend=c("EE", "EL","LL"),pch=1:3, border = NA)
}

# library(ggplot2)
# ggplot(df, aes(x = x, y = X.scaled, color = sex)) +
#   scale_x_continuous(name="Proximate cue", limits=range(df$x)) +
#   scale_y_continuous(name="Marine growth (scaled)", limits=range(df$y)) +
#   theme_bw() +
#   #geom_point() +
#   geom_point(df, aes(shape=genotype, color=sex)) +
#   stat_ellipse(aes(x=x, y=X.scaled,color=sex, group=sex),type = "norm") +
#   theme(legend.position='none')

#dev.off()












### THETAS & ETAS ####
#load("results/vgll3_scorff_EtasThetas_nimble.RData")


############### FIG 1


#mcmc <- do.call(rbind,samples)

## pool theta samples between the three chains
etas <- mcmc[,paste0("eta[",1:data$N,"]")]
thetas <- mcmc[,paste0("theta[",1:data$N,"]")]

etas.median <- apply(etas,2,median)
thetas.median <- apply(thetas,2,median)

#smp<-sample(1:dim(etas)[1],100,replace=FALSE)
smp <- ceiling(seq(1000,nrow(etas),10))
#smp.etas=smp.thetas=array(,dim=c(length(x),dim(etas)[2]))
#for (i in 1:length(x)){
smp.etas <- etas[smp,]
smp.thetas <- thetas[smp,]
#}


#### CORRELATIONS
cor.all=NULL
cor<-array(,dim=c(dim(etas)[1],3,2))
for (i in 1:dim(etas)[1]){
  cor.all[i] <-cor.test(etas[i,],thetas[i,])$estimate
  for (x in 1:2){
  for (gen in 1:3){
  tmp <- cbind(eta = etas[i,g==gen & sex==x],theta = thetas[i,g==gen & sex==x])
  cor[i,gen,x]<-cor.test(tmp[,1],tmp[,2])$estimate
  }}}
par(mfrow=c(1,1))
plot(NULL,xlim=c(0,5),ylim=c(-1,1),xlab="GENOTYPE",ylab="Correlation eta/theta", xaxt="n")
legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex)
axis(1, line=1, at=1:4,labels=c("All","EE", "EL", "LL"),tick = FALSE,cex.axis=1.5)
abline(h=0,lty=2)
segments(1, quantile(cor.all,probs=0.025),1, quantile(cor.all,probs=0.975),col=1)
segments(1, quantile(cor.all,probs=0.25),1, quantile(cor.all,probs=0.75),col=1,lwd=2)
points(1,quantile(cor.all,probs=0.5), pch=21, col=1,bg="white")

for (gen in (1:3)){
segments(gen-0.1+1, apply(cor[,,1], 2, quantile, probs=0.025)[gen],gen-0.1+1, apply(cor[,,1], 2, quantile, probs=0.975)[gen],col=col.sex[1])
segments(gen-0.1+1, apply(cor[,,1], 2, quantile, probs=0.25)[gen],gen-0.1+1, apply(cor[,,1], 2, quantile, probs=0.75)[gen],col=col.sex[1],lwd=2)
points(gen-0.1+1,apply(cor[,,1], 2, quantile, probs=0.5)[gen], pch=21, col=col.sex[1],bg="white")

segments(gen+0.1+1, apply(cor[,,2], 2, quantile, probs=0.025)[gen],gen+0.1+1, apply(cor[,,2], 2, quantile, probs=0.975)[gen],col=col.sex[2])
segments(gen+0.1+1, apply(cor[,,2], 2, quantile, probs=0.25)[gen],gen+0.1+1, apply(cor[,,2], 2, quantile, probs=0.75)[gen],col=col.sex[2],lwd=2)
points(gen+0.1+1,apply(cor[,,2], 2, quantile, probs=0.5)[gen], pch=21, col=col.sex[2],bg="white")
}


par(mfrow=c(1,1))




# smp<-sample(1:dim(etas)[1],1,replace=FALSE)
# 
# smp.etas <- smp.thetas <- NULL
# smp.etas <- etas[smp,]  
# smp.thetas <- thetas[smp,] 
# smp.cor <- cor[smp,]

par(mfrow=c(2,2))
for (x in 1:2){

smp=sample(1000:nrow(etas),1,replace=FALSE)  
plot(NULL,xlim=c(-7,7),ylim=c(-7,7),xlab="Proximate cue",ylab="Threshold", main=paste0("Male (iteration: ", smp,")"))
#legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex)
legend("topright", legend=c("EE", "EL", "LL"),pch=1:3, bty="n")
for (gen in 1:3){
  tmp <- cbind(eta = etas[smp,g==gen & sex==1],theta = thetas[smp,g==gen & sex==1])
  add_hd_ellipse(tmp, coverage = 0.95, border=NA, fill = "#63636370", lwd=3)
  points(etas[smp,g==gen & sex==1],thetas[smp,g==gen & sex==1], col=col.sex[1],pch=gen)
  #points(etas[smp,g==gen & sex==2],thetas[smp,g==gen & sex==2], col=col.sex[2],pch=gen)
}


plot(NULL,xlim=c(-5,5),ylim=c(-5,5),xlab="Proximate cue",ylab="Threshold",  main=paste0("Female (iteration: ", smp,")"))
smp=sample(1000:nrow(etas),1,replace=FALSE)
#legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex)
legend("topright", legend=c("EE", "EL", "LL"),pch=1:3, bty="n")
for (gen in 1:3){
  tmp <- cbind(eta = etas[smp,g==gen & sex==2],theta = thetas[smp,g==gen & sex==2])
  add_hd_ellipse(tmp, coverage = 0.95, border=NA, fill = "#63636370", lwd=3)
  points(etas[smp,g==gen & sex==2],thetas[smp,g==gen & sex==2], col=col.sex[2],pch=gen)
  #points(etas[smp,g==gen & sex==2],thetas[smp,g==gen & sex==2], col=col.sex[2],pch=gen)
}
}

# 
# #X=NULL
# #X <- data$Pronotum
# #
# par(mfcol=c(3,1))
# gene.name=c("EE","EL","LL")
# 
#   for (gene in 1:3){
#   #tmp <- subset(smp.etas[1,g==1], g==1)
#   plot(density(smp.etas[1,],adjust=2),bty="n",type='n', main=gene.name[gene],xlim=c(-10,10),ylim=c(0,10),xlab='',ylab='')
#   mtext("Density",side=2,line=2,at=.5,cex=1)
#   for (sexe in 1:2){
#   for (i in 1:length(smp)){
#     lines(density(smp.thetas[i,g==gene & sex==sexe]),col=colors[gene])
#     lines(density(smp.etas[i,g==gene  & sex==sexe]),col=t_col(col.sex[sexe],20))
#   }
#   }
#   }



############### FIG 2
par(mfcol=c(1,3))

ylim=c(-5,5)
xlim=c(-0.7,0.7)
smp=sample(1:nrow(etas),10,replace=FALSE)
smp.etas <- smp.thetas <- NULL
smp.etas <- etas[smp,]  
smp.thetas <- thetas[smp,] 

for (gene in 1:3){
  
  side=c(-1,1)
  #Make the plot
  if(gene==1){
    plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n")#,main=gene.name[gene])
    legend("topleft", legend=c("EE", "EL", "LL"),fill=colors, bty="n",border = colors, title ="Thresholds distributions")
    axis(1, line=-1, at=0,labels=c("EE"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
  } 
  if(gene==2){
    plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n",yaxt="n")#,main=gene.name[gene])
    legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex, title ="Cue distributions")
    axis(1, line=-1, at=0,labels=c("EL"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
    axis(1, line=2, at=0,labels=c("GENOTYPE"),tick = FALSE,cex.axis=2)#, "EL", "LL"))
  } 
  
  if(gene==3){
    plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n", yaxt="n")#,main=gene.name[gene])
    axis(1, line=-1, at=0,labels=c("LL"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
  }
  #axis(side = 1, at=0,labels = type)
  #axis(2, at=1:length(seq(40,140,10)),labels=seq(40,140,10), col.axis="red", las=2)
  
  text(-0.2,10,male,cex=2)
  text(0.2,10,female,cex=2)
  
  
  
  for (sexe in 1:2){
  #if(!isEmpty(smp.etas[sex==1 & gene==1])){
    for (i in 1:length(smp)){
    tmp = smp.etas[i,g==gene & sex==sexe]
    points(rnorm(length(tmp),0.5*side[sexe],0.025), tmp,pch=".",col=t_col(col.sex[sexe],10) )
    polygon(density(tmp)$y * side[sexe] * 2,density(tmp)$x , main="" , xlab="" , xaxt="n", col=NA , border=t_col(col.sex[sexe],50))
    }
    
    for (i in 1:length(smp)){
    tmp = smp.thetas[i,g==gene & sex==sexe]
    points(rnorm(length(tmp),0.5*side[sexe],0.025), tmp,pch=".",col=t_col(colors[gene],10) )
    polygon(density(tmp)$y * side[sexe] *0.2 ,density(tmp)$x , main="" , xlab="" , xaxt="n", col=NA , border=t_col(colors[gene],50))
    #segments(0,median(tmp),max(density(tmp)$y) * side[sexe],median(tmp),col="#99999950")
    }
  }
  }




############### FIG 3

# mu_alpha=NULL
# mu_alpha[1] <- 0.391995260 
# mu_alpha[2] <- 3.719397127
# 
# alpha=array(,dim=c(3,2))
# alpha[1, 1] <- -3.006811083
# alpha[2, 1] <- -2.940781596
# alpha[3, 1] <- 3.006811083
# alpha[1, 2] <- -2.750615029
# alpha[2, 2] <- -0.825334996
# alpha[3, 2] <-  2.750615029
# 
# sigma2_theta=array(,dim=c(3,2))
# sigma2_theta[1, 1]  <- 1.415155727
# sigma2_theta[2, 1]  <-  1.415155727
# sigma2_theta[3, 1]  <-  1.415155727 
# sigma2_theta[1, 2]  <-  0.327978080
# sigma2_theta[2, 2]  <-  0.327978080
# sigma2_theta[3, 2]  <-  0.327978080
# 
# sigma2_eta=NULL
# sigma2_eta[1] <-  9.553079998
# sigma2_eta[2] <- 9.553079998
# 
# #etas <- rnorm(100,0, sqrt(1+sigma2_eta[sexe]))
# 
# n=10000
# thetas=etas=array(,dim=c(n,3,2))
# 
# par(mfcol=c(1,3))
# #par(mfcol=c(1,3),
# #oma=c(1,1,0,0),
# #mar=c(1,1,1,0),
# #mgp=c(3,1,0)
# #xpd=NA
# #)
# ylim=c(-15,17)
# xlim=c(-0.7,0.7)
# 
# for (gene in 1:3){
#   
#   side=c(-1,1)
#   #Make the plot
#   if(gene==1){
#     plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n")#,main=gene.name[gene])
#     legend("topleft", legend=c("EE", "EL", "LL"),fill=colors, bty="n",border = colors, title ="Thresholds distributions")
#     axis(1, line=0, at=0,labels=c("EE"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
#      } 
#   if(gene==2){
#     plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n",yaxt="n")#,main=gene.name[gene])
#     legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex, title ="Cue distributions")
#     axis(1, line=0, at=0,labels=c("EL"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
#     axis(1, line=2, at=0,labels=c("GENOTYPE"),tick = FALSE,cex.axis=2)#, "EL", "LL"))
#     } 
#     
#   if(gene==3){
#     plot(NULL, ylim=ylim,xlim=xlim,bty="n",xlab="",ylab="", xaxt="n", yaxt="n")#,main=gene.name[gene])
#     axis(1, line=0, at=0,labels=c("LL"),tick = FALSE,cex.axis=1.5)#, "EL", "LL"))
#     }
#   #axis(side = 1, at=0,labels = type)
#   #axis(2, at=1:length(seq(40,140,10)),labels=seq(40,140,10), col.axis="red", las=2)
#   
#   text(-0.2,-15,male,cex=2)
#   text(0.2,-15,female,cex=2)
#   
#   
#   for (sexe in 1:2){
#     
#     thetas[,gene,sexe] <- rnorm(n,mu_alpha[sexe] + alpha[gene, sexe], sqrt(sigma2_theta[gene, sexe]))
#     etas[,gene,sexe] <- rnorm(n,0, sqrt(1+sigma2_eta[sexe]))
#     
#     #if(!isEmpty(smp.etas[sex==1 & gene==1])){
#     tmp = etas[ ,gene, sexe]
#     tmp2= etas.median[data$g==gene & data$sex==sexe]
#     tmp3= thetas.median[data$g==gene & data$sex==sexe]
#     
#     points(rnorm(length(tmp),0.5*side[sexe],0.025), tmp,pch=".",col=t_col(col.sex[sexe],10) )
#     #points(rnorm(length(tmp2),0.5*side[sexe],0.025),tmp2,pch=16,col=t_col(col.sex[sexe],20) )
#     #points(rnorm(length(tmp3),0.5*side[sexe],0.025),tmp3,pch=16,col=t_col(colors[gene],20) )
#     
#     
#     polygon(density(tmp)$y * side[sexe],density(tmp)$x , main="" , xlab="",ylim=c(-10,10) , xaxt="n", col=t_col(col.sex[sexe],20) , border=NA)
#     segments(0,median(tmp),max(density(tmp)$y) * side[sexe],median(tmp),col="#99999950")
# 
#     tmp = thetas[ ,gene, sexe]
#     points(rnorm(length(tmp),0.5*side[sexe],0.025), tmp,pch=".",col=t_col(colors[gene],10) )
#     
#     polygon(density(tmp)$y * side[sexe],density(tmp)$x , main="" , xlab="",ylim=c(-10,10) , xaxt="n", col=t_col(colors[gene],50) , border=NA)
#    segments(0,median(tmp),max(density(tmp)$y) * side[sexe],median(tmp),col="#99999950")
#     
#   }
# }











#load("~/Documents/RESEARCH/PROJECTS/LETM/VGLL3/results/vgll3_scorff_jags.RData")

#thetas <- samples$BUGSoutput$sims.list$theta
#etas <- (samples$BUGSoutput$sims.list$eta)

#thetas_male <- thetas[,sex==1]
#thetas_female <- thetas[,sex==2]


#thetas_medians=array(,dim=c(nrow(thetas), 33,2))
thetas_medians=array(,dim=c(5, nyears,2));rownames(thetas_medians)<-c("2.5%","25%","50%","75%","97.5%")
etas_medians=array(,dim=c(5, nyears,2));rownames(etas_medians)<-c("2.5%","25%","50%","75%","97.5%")
etas_medians_all=array(,dim=c(5, nyears));rownames(etas_medians_all)<-c("2.5%","25%","50%","75%","97.5%")

  # for (y in 1987:2016){
  #   #thetas_medians_all[,y-1985,s] <- quantile(apply(thetas[,sex==s & year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
  #   etas_medians_all[,y-1985] <- quantile(apply(etas[,year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
  #   for (s in 1:2){
  #   thetas_medians[,y-1985,s] <- quantile(apply(thetas[,sex==s & year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
  #   #etas_medians[,y-1985,s] <- quantile(apply(etas[,sex==s & year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
  # }              
  # }

years <- sort(unique(data$year))
nyears <- length(unique(data$year))
etas_means_all=array(,dim=c(5, nyears));rownames(etas_means_all)<-c("2.5%","25%","50%","75%","97.5%")
colnames(etas_means_all)<-years

thetas_means_male=array(,dim=c(5, nyears))
rownames(thetas_means_male)<-c("2.5%","25%","50%","75%","97.5%")
colnames(thetas_means_male)<-years
thetas_means_female=etas_means_female=etas_means_male=thetas_means_male

j=0
for (y in years){
  j=j+1
  #print(j)
  tmp=NULL
  tmp <-etas[,which(year==y)]
  if(is.null(ncol(tmp))==FALSE) {
    etas_means_all[,j] <- quantile(apply(tmp,1,mean),probs=c(0.025,0.25, 0.5,0.75, 0.975))
 
tmp_male <- thetas[,which(sex==1 & year==y)]
thetas_means_male[,j] <- quantile(apply(tmp_male,1,mean),probs=c(0.025,0.25, 0.5,0.75, 0.975))
tmp_male <- etas[,which(sex==1 & year==y)]
etas_means_male[,j] <- quantile(apply(tmp_male,1,mean),probs=c(0.025,0.25, 0.5,0.75, 0.975))

tmp_female <- thetas[,which(sex==2 & year==y)]
thetas_means_female[,j] <- quantile(apply(tmp_female,1,mean),probs=c(0.025,0.25, 0.5,0.75, 0.975))
tmp_female <- etas[,which(sex==2 & year==y)]
etas_means_female[,j] <- quantile(apply(tmp_female,1,mean),probs=c(0.025,0.25, 0.5,0.75, 0.975))
  }
}



# thetas_50 <- samples$BUGSoutput$median$theta
# etas_50 <- samples$BUGSoutput$median$eta
# thetas_medians_ind=etas_medians_ind=array(,dim=c(data$N,32,2))
# for (s in 1:2){
#   for (y in 1987:2016){
#     tmp <- thetas_50[sex==s & year==y]
#     thetas_medians_ind[1:length(tmp),y-1985,s] <- tmp
#     tmp=NULL
#     tmp <- etas_50[sex==s & year==y]
#     etas_medians_ind[1:length(tmp),y-1985,s] <- tmp
#   }
# }

thetas_means_ind=etas_means_ind=array(,dim=c(data$N,nyears,2))
for (s in 1:2){
  for (y in years){
    tmp <- thetas_means[sex==s & year==y]
    thetas_means_ind[1:length(tmp),y-min(data$year)-1,s] <- tmp
    tmp=NULL
    tmp <- etas_means[sex==s & year==y]
    etas_means_ind[1:length(tmp),y-min(data$year)-1,s] <- tmp
  }
}



#par(mfcol=c(1,2))

# #Xc <- (X-mean(X))/sd(X)
# X.scaled=NULL
# for (i in 1: length(X)){
#   X.scaled[i] <- (X[i]-mu_X[data$g[i],data$sex[i]])/sd_X[data$g[i],data$sex[i]]
# }
# Xc <- X.scaled
# 
# #tmp <- aggregate(X~year+sex, FUN = mean)
# tmp <- aggregate(Xc~year+sex, FUN = 'quantile', probs=c(0.025,0.25, 0.5,0.75, 0.975))#;colnames(tmp)<-c("year","sex","2.5%","50%","97.5%")
# 
# range_years <- 1986:2017
# 
# plot(NULL, xlim=c(1985,2017),ylim=c(-3,3), ylab="Growth (scaled)")
# points(range_years-0.1, tmp[,"Xc"][tmp$sex==1,"50%"],pch=16,col=1)
# segments(range_years-.1,tmp[,"Xc"][tmp$sex==1,"2.5%"], range_years-.1,tmp[,"Xc"][tmp$sex==1,"97.5%"],col=1)
# segments(range_years-.1,tmp[,"Xc"][tmp$sex==1,"25%"], range_years-.1,tmp[,"Xc"][tmp$sex==1,"75%"],col=1,lwd=2)
# points(range_years+0.1, tmp[,"Xc"][tmp$sex==2,"50%"],pch=16,col=2)
# segments(range_years+.1,tmp[,"Xc"][tmp$sex==2,"2.5%"], range_years+.1,tmp[,"Xc"][tmp$sex==2,"97.5%"],col=2)
# segments(range_years+.1,tmp[,"Xc"][tmp$sex==2,"25%"], range_years+.1,tmp[,"Xc"][tmp$sex==2,"75%"],col=2,lwd=2)
# 
# 
# y <- tmp[,"Xc"][tmp$sex==1,"50%"]
# x <- range_years
# mod1=loess(y~x,span=0.5)
# xfit=seq(from=min(x),to=max(x),length.out=30)
# yfit1=predict(mod1,newdata=xfit)
# points(xfit,yfit1,type="l",lwd=2,col=1)
# 
# y <- tmp[,"Xc"][tmp$sex==2,"50%"]
# x <- range_years
# mod1=loess(y~x,span=0.5)
# xfit=seq(from=min(x),to=max(x),length.out=30)
# yfit1=predict(mod1,newdata=xfit)
# points(xfit,yfit1,type="l",lwd=2,col=2)



#par(mfcol=c(2,2))
# Define the layout matrix
# layout_matrix <- matrix(c(1, 2, 3,3), nrow = 2, ncol = 2, byrow = TRUE)
# # Set up the layout
# layout(layout_matrix)
par(mfcol=c(1,2))
range_years <- as.numeric(colnames(etas_means_all))
plot(NULL, xlim=range(range_years),ylim=c(-1,1), ylab="Proximate cue",xlab="")
#rect(xleft = 1985, xright = 2005, ybottom = -4, ytop = 4,border = "lightgrey", col = "lightgrey")
#rect(xleft = 2005, xright = 2018, ybottom = -4, ytop = 4,border = "darkgrey", col = "darkgrey")

segments(range_years-.1,etas_means_all["2.5%",], range_years-.1,etas_means_all["97.5%",] ,col=1)
segments(range_years-.1,etas_means_all["25%",], range_years-.1,etas_means_all["75%",] ,col=1,lwd=2)
points(range_years-.1,etas_means_all["50%",],col=1,pch=16)
# for (y in 1:32){
#   for (i in 1:nrow(etas_medians_ind[,,1])){
#  #   points(1985+y, etas_medians_ind[i,y,1],col=1,pch=3)
#   }
#   }

y <- etas_means_all["50%",]
x <- range_years
mod1=loess(y~x,span=0.5)
xfit=seq(from=min(x),to=max(x),length.out=30)
yfit1=predict(mod1,newdata=xfit)
points(xfit,yfit1,type="l",lwd=2,col=1)

# segments(range_years-.1,etas_medians["2.5%",,1], range_years-.1,etas_medians["97.5%",,1] ,col=1)
# segments(range_years-.1,etas_medians["25%",,1], range_years-.1,etas_medians["75%",,1] ,col=1,lwd=2)
# points(range_years-.1,etas_medians["50%",,1],col=1,pch=16)
# # for (y in 1:32){
# #   for (i in 1:nrow(etas_medians_ind[,,1])){
# #  #   points(1985+y, etas_medians_ind[i,y,1],col=1,pch=3)
# #   }
# #   }
# 
# y <- etas_medians["50%",,1]
# x <- range_years
# mod1=loess(y~x,span=0.5)
# xfit=seq(from=min(x),to=max(x),length.out=30)
# yfit1=predict(mod1,newdata=xfit)
# points(xfit,yfit1,type="l",lwd=2,col=1)
# 
# 
# segments(range_years+.1,etas_medians["2.5%",,2], range_years+.1,etas_medians["97.5%",,2] ,col=2)
# segments(range_years+.1,etas_medians["25%",,2], range_years+.1,etas_medians["75%",,2] ,col=2,lwd=2)
# points(range_years+.1,etas_medians["50%",,2],col=2,pch=16)
# # for (y in 1:32){
# #   for (i in 1:nrow(etas_medians_ind[,,2])){
# # #    points(1985+y, etas_medians_ind[i,y,2],col=2,pch=3)
# #   }
# # }
# y <- etas_medians["50%",,2]
# x <- range_years
# mod1=loess(y~x,span=0.5)
# xfit=seq(from=min(x),to=max(x),length.out=30)
# yfit1=predict(mod1,newdata=xfit)
# points(xfit,yfit1,type="l",lwd=2,col=2)





#par(mfcol=c(1,1))
#tmp <- aggregate(thetas_50~year+sex, FUN = mean)
range_years <- as.numeric(colnames(thetas_means_male))
plot(NULL, xlim=range(range_years),ylim=c(-6,4), ylab="Tresholds",xlab="")
# rect(xleft = 1986, xright = 2005, ybottom = -5, ytop = 5,border = "lightgrey", col = "lightgrey")
# rect(xleft = 2005, xright = 2016, ybottom = -5, ytop = 5,border = "darkgrey", col = "darkgrey")

# for (y in 1986:2017){
# points(y-.1,thetas_50[year==y & ],col=1,pch=16)
# }

segments(range_years-.1,thetas_means_male["2.5%",], range_years-.1,thetas_means_male["97.5%",] ,col=1)
segments(range_years-.1,thetas_means_male["25%",], range_years-.1,thetas_means_male["75%",] ,col=1,lwd=2)
points(range_years-.1,thetas_means_male["50%",],col=1,pch=16)
# for (y in 1:32){
#   for (i in 1:nrow(thetas_medians_ind[,,1])){
# #    points(1985+y, thetas_medians_ind[i,y,1],col=1,pch=3)
#   }
# }

y <- thetas_means_male["50%",]
x <- range_years
mod1=loess(y~x,span=0.5)
xfit=seq(from=min(x),to=max(x),length.out=30)
yfit1=predict(mod1,newdata=xfit)
points(xfit,yfit1,type="l",lwd=2,col=1)


segments(range_years-.1,thetas_means_female["2.5%",], range_years-.1,thetas_means_female["97.5%",] ,col=2)
segments(range_years-.1,thetas_means_female["25%",], range_years-.1,thetas_means_female["75%",] ,col=2,lwd=2)
points(range_years-.1,thetas_means_female["50%",],col=2,pch=16)
for (y in 1:nyears){
  for (i in 1:nrow(thetas_means_ind[,,2])){
 #   points(1985+y, thetas_medians_ind[i,y,2],col=2,pch=3)
  }
}

y <- thetas_means_female["50%",]
x <- range_years
mod1=loess(y~x,span=0.5)
xfit=seq(from=min(x),to=max(x),length.out=30)
yfit1=predict(mod1,newdata=xfit)
points(xfit,yfit1,type="l",lwd=2,col=2)

legend("topright",legend=c(male, female),text.col =c(1,2), bty="n")










#linearMod <- lm(y ~ x)  # build linear regression model on full data
#print(linearMod)

par(mfcol=c(1,3))
plot(thetas_means_male["50%",],thetas_means_female["50%",],pch="", xlab="MALE",ylab="FEMALE", main="Tresholds (medians)")
text(thetas_means_male["50%",],thetas_means_female["50%",], min(data$year):2016)
# Compute correlation
cor_result <- cor.test(thetas_means_male["50%",], thetas_means_female["50%",])

# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)

# Add correlation text to the plot (top-left corner)
text(x=min(thetas_means_male["50%",]), 
     y=max(thetas_means_female["50%",])*0.95, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="tomato")  # 'pos=4' places text to the right




plot(etas_means_all["50%",],thetas_means_male["50%",],pch="", xlab="eta",ylab="theta", main="MALE")
text(etas_means_all["50%",],thetas_means_male["50%",], 1986:2016)
# Compute correlation
cor_result <- cor.test(etas_means_all["50%",], thetas_means_male["50%",])
# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)
# Add correlation text to the plot (top-left corner)
text(x=min(etas_means_all["50%",]), 
     y=max(thetas_means_male["50%",])*0.9, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="tomato")  # 'pos=4' places text to the right


plot(etas_means_all["50%",],thetas_means_female["50%",],pch="", xlab="eta",ylab="theta", main="FEMALE")
text(etas_means_all["50%",],thetas_means_female["50%",], 1986:2016)
# Compute correlation
cor_result <- cor.test(etas_means_all["50%",], thetas_means_female["50%",])
# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)
# Add correlation text to the plot (top-left corner)
text(x=min(etas_means_all["50%",]), 
     y=max(thetas_means_female["50%",])*0.95, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="tomato")  # 'pos=4' places text to the right





par(mfcol=c(1,1))
plot(etas_means_male["50%",1:nyears],thetas_means_male["50%",1:nyears], xlim=c(-1,1), ylim=c(-3.5,2), pch="",main="Male (black) / Female (red)", xlab="Proximate cues",ylab="Thresholds")
points(etas_means_male["50%",1:nyears],thetas_means_male["50%",1:nyears], pch="")
text(etas_means_male["50%",1:nyears],thetas_means_male["50%",1:nyears],1986:2016,cex=0.75)
# Compute correlation
cor_result <- cor.test(etas_means_male["50%",], thetas_means_male["50%",])
# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)
# Add correlation text to the plot (top-left corner)
text(x=-1, 
     y=max(thetas_means_male["50%",])*0.9, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="black")  # 'pos=4' places text to the right

#plot(etas_means_all["50%",1:nyears],thetas_means_female["50%",1:nyears], xlim=c(-1,1), ylim=c(-0.5,2), pch="",main="Female", xlab="Proximate cues",ylab="Thresholds")
points(etas_means_female["50%",1:nyears],thetas_means_female["50%",1:nyears], pch="")
text(etas_means_female["50%",1:nyears],thetas_means_female["50%",1:nyears],1986:2016,cex=0.75,col="tomato")
# Compute correlation
cor_result <- cor.test(etas_means_female["50%",], thetas_means_female["50%",])
# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)
# Add correlation text to the plot (top-left corner)
text(x=-1, 
     y=max(thetas_means_female["50%",])*0.9, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="tomato")  # 'pos=4' places text to the right
legend("topright",legend=c(male, female),text.col =c(1,2), bty="n")





#### ALLELE FREQUENCIES
#load("data/data_vgll3.rdata")
#years<-sort(unique(df$t))
#nyears<-length(years)
par(mfcol=c(1,1))
tmp <- aggregate(X~g+sex+year, FUN=length)
n <- xtabs( X ~ g+sex+year, tmp) # convert dataframe to matrix
## Frequence of allele E
p=q=array(,dim=c(2,nyears))
for (s in 1:2){    # sex
  for (t in 1:nyears){
    p[s,t] <- sum((2*n[1,s,t])+(1*n[2,s,t])+(0*n[3,s,t]))/(2*sum(n[1:3,s,t])) # population frequence allelic for E
    q[s,t] <- 1 - p[s,t]
  }}
colnames(p)<-years;rownames(p)<-c("Male","Female")

## Load proportion of males by age at sea (1SW vs MSW, and by year)
load("data/p_male_Scorff.Rdata")
p_male_1SW <- mean(p_male_1SW[,"50%"])
p_female_1SW <- 1-p_male_1SW
p_male_MSW <- mean(p_male_MSW[,"50%"])
p_female_MSW <- 1-p_male_MSW


# loading coda
#load(paste('results/Results_',stade,"_",year,'.RData',sep=""))
#fit.mcmc <- as.mcmc(fit$sims.matrix) # using bugs
#n1SW <- as.matrix(fit.mcmc[,paste0("n_1SW[",1:data$Y,"]")])
#ntot <- as.matrix(fit.mcmc[,paste0("n_tot[",1:data$Y,"]")])
p_1SW <- 0.86 #median(n1SW/ntot)

P <- array(,dim=c(2,nyears))
P[1,] <- q[1,]*p_1SW*p_male_1SW + q[1,]*(1-p_1SW)*p_male_MSW
P[2,] <- q[2,]*p_1SW*p_female_1SW + q[2,]*(1-p_1SW)*p_female_MSW

col.sex=c(1,2)
par(mfcol=c(1,1))
plot(NULL,xlim=c(min(years),max(years)),ylim=c(0,.3),xlab="",ylab="Allele L frequencies")
for (s in 1:2){
  y <- P[s,][-1]
  x <- sort(years)[-1]
  mod1=loess(y~x,span=0.9)
  xfit=seq(from=min(x),to=max(x),length.out=nyears-1)
  yfit1=predict(mod1,newdata=xfit)
  points(xfit,yfit1,type="l",lwd=2,col=col.sex[s])
  points(years,P[s,],col=col.sex[s],pch=16)
}
legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex, title ="Allele L frequencies")



#Expected Heterozygosity 
He <- array(,dim=c(2,nyears))
He[1,] <- 1- ((q[1,]*p_1SW*p_male_1SW + q[1,]*(1-p_1SW)*p_male_MSW)^2 + (p[1,]*p_1SW*p_male_1SW + p[1,]*(1-p_1SW)*p_male_MSW)^2)
#P[2,] <- q[2,]*p_1SW*p_female_1SW + q[2,]*(1-p_1SW)*p_female_MSW
He[2,] <- 1- ((q[2,]*p_1SW*p_female_1SW + q[2,]*(1-p_1SW)*p_female_MSW)^2 + (p[2,]*p_1SW*p_female_1SW + p[2,]*(1-p_1SW)*p_female_MSW)^2)

#Observed Heterozygosity 
he=array(,dim=c(2,nyears))
for (s in 1:2){    # sex
  for (t in 1:nyears){
    he[s,t] <- sum(n[2,s,t])/(sum(n[1:3,s,t])) # population frequence allelic for E
  }}
he[1,] <- he[1,] *p_1SW*p_male_1SW + he[1,] *(1-p_1SW)*p_male_MSW
he[2,] <- he[1,] *p_1SW*p_female_1SW + he[2,] *(1-p_1SW)*p_female_MSW

col.sex=c(1,2)
par(mfcol=c(1,1))
plot(NULL,xlim=c(min(years),max(years)),ylim=c(0,1),xlab="",ylab="Heterozozygosity")
for (s in 1:2){
  y <- He[s,][-1]
  x <- sort(years)[-1]
  mod1=loess(y~x,span=0.9)
  xfit=seq(from=min(x),to=max(x),length.out=nyears-1)
  yfit1=predict(mod1,newdata=xfit)
  points(xfit,yfit1,type="l",lwd=2,col=col.sex[s])
  points(years,He[s,],col=col.sex[s],pch=16)
}

for (s in 1:2){
  y <- he[s,][-1]
  x <- sort(years)[-1]
  mod1=loess(y~x,span=0.9)
  xfit=seq(from=min(x),to=max(x),length.out=nyears-1)
  yfit1=predict(mod1,newdata=xfit)
  points(xfit,yfit1,type="l",lwd=2,col=col.sex[s], lty=2)
  points(years,he[s,],col=col.sex[s],pch=17)
}
legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex, title ="")







par(mfcol=c(1,1))
plot(P[1,1:nyears],thetas_means_male["50%",], xlim=c(0, .3), ylim=c(-3.5,2), pch="",main="", xlab="Freq allele L",ylab="Thresholds")
points(P[1,1:nyears],thetas_means_male["50%",], pch="")
text(P[1,1:nyears],thetas_means_male["50%",],min(data$year):max(data$year),cex=0.75)
# Compute correlation
cor_result <- cor.test(P[1,1:nyears],thetas_means_male["50%",])
# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)
# Add correlation text to the plot (top-left corner)
text(x=0, 
     y=max(thetas_means_male["50%",])*0.95, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="black") 

#plot(P[2,1:nyears],thetas_means_female["50%",], xlim=c(0, .3), ylim=c(-3,2), pch="", main="Female", xlab="Freq allele L",ylab="Thresholds")
points(P[2,1:nyears],thetas_means_female["50%",], pch="")
text(P[2,1:nyears],thetas_means_female["50%",],min(data$year):max(data$year),cex=0.75,col="tomato")
# Compute correlation
cor_result <- cor.test(P[2,1:nyears],thetas_means_female["50%",])
# Extract correlation estimate (r) and p-value
cor_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)
# Add correlation text to the plot (top-left corner)
text(x=0, 
     y=max(thetas_means_female["50%",])*0.95, 
     labels=paste0("r = ", cor_value, "\np = ", p_value), 
     pos=4, col="tomato") 
legend("topright",legend=c(male, female),text.col =c(1,2), bty="n")





# 
# par(mfcol=c(1,1))
# 
# thetas <- samples$BUGSoutput$sims.list$theta
# etas <- samples$BUGSoutput$sims.list$eta
# 
# 
# thetas.stat <- apply(thetas,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# etas.stat <- apply(etas,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# 
# smp=sample(1:data$N,50,replace=FALSE)
# 
# plot(NULL, xlim=c(0,50), ylim=c(-10,10),xlab="",ylab="Individuals")
# for (i in 1:50){
#   segments(i-.1,thetas.stat["2.5%",smp[i]], i-.1,thetas.stat["97.5%",smp[i]] ,col=1)
#   segments(i-.1,thetas.stat["25%",smp[i]], i-.1,thetas.stat["75%",smp[i]] ,col=1,lwd=2)
#   points(i-0.1, thetas.stat[3,smp[i]], pch=21,bg="white",col=1)
#   
#   segments(i+.1,etas.stat["2.5%",smp[i]], i+.1,etas.stat["97.5%",smp[i]] ,col=2)
#   segments(i+.1,etas.stat["25%",smp[i]], i+.1,etas.stat["75%",smp[i]] ,col=2,lwd=2)
#   points(i+.1, etas.stat[3,smp[i]], pch=21,bg="white",col=2)
#   
#   #points(i+.1, data$X.scaled[smp[i]], pch="X",bg="white",col=2)
#   points(i+.1, X.scaled[smp[i]], pch="X",bg="white",col=2)
#   text(i,-10,labels=data$Y[smp[i]],cex=0.75)  
# }
# legend("topleft",legend=c("Threshold","Prox. Cue","Growth scaled"),pch=c(1,1,"X"),col=c(1,2,2) ,bty="n")


dev.off()


