rm(list=ls())   # Clear memory

#____________________DATA______________________#
source("code/DATA_VGLL3_Scorff.R")
attach(dataToJags)


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


nimble=FALSE
jags=TRUE


# ---- DISTRIBUTIONS -----
if(nimble){
  load("results/vgll3_scorff_EtaSameSex_nimble.RData")
  stats <- MCMCpstr(samples, func = function(x) quantile(x, probs = c(0.025,0.5 ,0.975)))
  mcmc <- do.call(rbind,samples)
  
  ## pool theta samples between the three chains
  alpha1m <- mcmc[,"alpha[1, 1]"]
  #alpha1m <- c(samples$chain1[,"alpha[1, 1]"], samples$chain2[,"alpha[1, 1]"],samples$chain3[,"alpha[1, 1]"])
  alpha2m <- c(samples$chain1[,"alpha[2, 1]"], samples$chain2[,"alpha[2, 1]"],samples$chain3[,"alpha[2, 1]"])
  alpha3m <- c(samples$chain1[,"alpha[3, 1]"], samples$chain2[,"alpha[3, 1]"],samples$chain3[,"alpha[3, 1]"])
  
  
  alpha1f <- c(samples$chain1[,"alpha[1, 2]"], samples$chain2[,"alpha[1, 2]"],samples$chain3[,"alpha[1, 2]"])
  alpha2f <- c(samples$chain1[,"alpha[2, 2]"], samples$chain2[,"alpha[2, 2]"],samples$chain3[,"alpha[2, 2]"])
  alpha3f <- c(samples$chain1[,"alpha[3, 2]"], samples$chain2[,"alpha[3, 2]"],samples$chain3[,"alpha[3, 2]"])
  
  mu_alphaMale <- c(samples$chain1[,"mu_alpha[1]"], samples$chain2[,"mu_alpha[1]"],samples$chain3[,"mu_alpha[1]"])
  mu_alphaFemale <- c(samples$chain1[,"mu_alpha[2]"], samples$chain2[,"mu_alpha[2]"],samples$chain3[,"mu_alpha[2]"])
  
  alphas <- stats$alpha
  mu_alphas <- stats$mu_alpha
  
  sigma2_eta <- stats$sigma2_eta
  sigma2_theta <- stats$sigma2_theta
  sigma2_THETA <- stats$sigma2_THETA
}

if(jags){
  load("results/RESULTS_vgll3_scorff.RData")
  mcmc <- samples$BUGSoutput$sims.matrix
  #mu_alphaMale <- mcmc[,"mu_theta[1,1]"]#samples$BUGSoutput$sims.list$mu_theta[,1]
  #mu_alphaFemale <- mcmc[,"mu_theta[1,2]"]#samples$BUGSoutput$sims.list$mu_theta[,2]
  mu_alphaMale <- samples$BUGSoutput$sims.list$mu_theta[,1]
  mu_alphaFemale <- samples$BUGSoutput$sims.list$mu_theta[,2]
  alpha1m <- mcmc[,"alpha[1,1]"]
  alpha2m <- mcmc[,"alpha[2,1]"]
  alpha3m <- mcmc[,"alpha[3,1]"]
  
  alpha1f <- mcmc[,"alpha[1,2]"]
  alpha2f <- mcmc[,"alpha[2,2]"]
  alpha3f <- mcmc[,"alpha[3,2]"]

  theta <- array(,dim=c(nrow(mu_alphaFemale),3,2))
  theta[,1,1] <- mu_alphaMale + alpha1m
  theta[,2,1] <- mu_alphaMale + alpha2m
  theta[,3,1] <- mu_alphaMale + alpha3m
  
  theta[,1,2] <- mu_alphaFemale + alpha1f
  theta[,2,2] <- mu_alphaFemale + alpha2f
  theta[,3,2] <- mu_alphaFemale + alpha3f
  
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
}



pdf("results/FIGURES_VGLL3_Scorff.pdf")

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
ylim=c(1,4.5)
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




### ETA
### ## TODO




par(mfrow=c(2,2))
#boxplot(X~g + sex, data=new.df)
plot(NULL, xlim=c(0,4),ylim=c(2,4),xlab="Genotype", ylab="Growth",xaxt="n")
axis(1, at=1:3,labels=c("EE", "EL", "LL"))
legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")

for (j in 1:3){
  # MALE
  segments(j-0.1, quantile(new.df$X[g==j & sex==1], probs = 0.025),j-0.1, quantile(new.df$X[g==j & sex==1], probs = 0.975),col=colors[j])
  segments(j-0.1, quantile(new.df$X[g==j & sex==1], probs = 0.25),j-0.1, quantile(new.df$X[g==j & sex==1], probs = 0.75),col=colors[j],lwd=2)
  points(j-0.1,median(new.df$X[g==j & sex==1]), pch=24, col=colors[j],bg="white")
  text(j-0.1,quantile(new.df$X[g==j & sex==1], probs = 0.995),male,col=colors[j])

  # FEMALE
  segments(j+0.1, quantile(new.df$X[g==j & sex==2], probs = 0.025),j+0.1, quantile(new.df$X[g==j & sex==2], probs = 0.975),col=colors[j])
  segments(j+0.1, quantile(new.df$X[g==j & sex==2], probs = 0.25),j+0.1, quantile(new.df$X[g==j & sex==2], probs = 0.75),col=colors[j],lwd=2)
  points(j+0.1,median(new.df$X[g==j & sex==2]), pch=21, col=colors[j],bg="white")
  text(j+0.1,quantile(new.df$X[g==j & sex==2], probs = 0.995),female,col=colors[j])
}



#par(mfrow=c(2,2))
plot(NULL, xlim=c(0,4),ylim=c(-3,3),xlab="Genotype", ylab="Means of maturation thresholds",xaxt="n")
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
#par(mfrow=c(1,2))
X.sim = seq(1,5,0.1)
X.pred <- (X.sim-mean(X))/sd(X) #scale(new.df$X, center = TRUE, scale = TRUE))

X.obs <- seq(min(X), max(X), 0.1)
X.obsc <- (X.obs-mean(X))/sd(X)

p.pred=array(,dim=c(length(X.pred), 3, 2))
p.obs=array(,dim=c(length(X.obs), 3, 2))
#mu_alphas <- samples$BUGSoutput$median$mu_alpha
#alphas <- samples$BUGSoutput$median$alpha
#sigma2_eta <- samples$BUGSoutput$median$sigma2_eta
#sigma2_THETA <- samples$BUGSoutput$median$sigma2_THETA

#thetas <- mu_alphas + alphas

#p[,1,1] <- pnorm((X.pred - alphas[1,1])/ (sigma_eta + sqrt(sigma2_THETA[1,1])))

#for (i in 1:length(X.pred)){

mu_X <- samples$BUGSoutput$median$mu_X
sigma_X <- sqrt(samples$BUGSoutput$median$sigma2_X)
sigma2_eta <- samples$BUGSoutput$median$sigma2_eta



for (i in 1:3){
  
  
  
  #z[i]<-(((X.scaled) / sqrt(sigma2_eta[g[i],sex[i]] + 1)) - theta[i]) / sqrt(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]] + 1))
  #z[i]<-(((X[i]-mu_X[g[i],sex[i]]) / sqrt(sigma2_eta[g[i],sex[i]] + sigma2_X[g[i],sex[i]])) - theta[i]) / sqrt(sigma2_eta[g[i],sex[i]]/(sigma2_eta[g[i],sex[i]] + sigma2_X[g[i],sex[i]]))
  
  X.scaled <- (X.sim-mu_X[i,1])/sigma_X[i,1]
  theta_med <- (apply(theta[,,1], 2, quantile, probs=0.5))
  z <- (((X.scaled) / sqrt(sigma2_eta[i,1] + 1)) - theta_med[i]) / sqrt(sigma2_eta[i,1]/(sigma2_eta[i,1] + 1))
  tmp <- pnorm(z)
  #tmp <- pnorm( (((X.sim - mu_X[i,1])/sqrt(sigma2_eta[i,1] + sigma2_X[i,1])) - (apply(theta[,,1], 2, quantile, probs=0.5)[i])) / sqrt(sigma2_eta[i,1] / (sigma2_eta[i,1] + sigma2_X[i,1])) )
  p.pred[, i, 1] <- as.vector(tmp)
  
  X.scaled <- (X.obs-mu_X[i,1])/sigma_X[i,1]
  theta_med <- (apply(theta[,,1], 2, quantile, probs=0.5))
  z <- (((X.scaled) / sqrt(sigma2_eta[i,1] + 1)) - theta_med[i]) / sqrt(sigma2_eta[i,1]/(sigma2_eta[i,1] + 1))
  tmp <- pnorm(z)
  #tmp <- pnorm( (((X.obs - mu_X[i,1])/sqrt(sigma2_eta[i,1] + sigma2_X[i,1])) - (apply(theta[,,1], 2, quantile, probs=0.5)[i])) / sqrt(sigma2_eta[i,1] / (sigma2_eta[i,1] + sigma2_X[i,1])) )
  p.obs[, i, 1] <- as.vector(tmp)
  
  
  X.scaled <- (X.sim-mu_X[i,2])/sigma_X[i,2]
  theta_med <- (apply(theta[,,2], 2, quantile, probs=0.5))
  z <- (((X.scaled) / sqrt(sigma2_eta[i,2] + 1)) - theta_med[i]) / sqrt(sigma2_eta[i,2]/(sigma2_eta[i,2] + 1))
  tmp <- pnorm(z)
  #tmp <- pnorm( (((X.sim - mu_X[i,2])/sqrt(sigma2_eta[i,2] + sigma2_X[i,2])) - (apply(theta[,,2], 2, quantile, probs=0.5)[i])) / sqrt(sigma2_eta[i,2] / (sigma2_eta[i,2] + sigma2_X[i,2])) )
  p.pred[, i, 2] <- as.vector(tmp)
  
  X.scaled <- (X.obs-mu_X[i,2])/sigma_X[i,2]
  theta_med <- (apply(theta[,,2], 2, quantile, probs=0.5))
  z <- (((X.scaled) / sqrt(sigma2_eta[i,2] + 1)) - theta_med[i]) / sqrt(sigma2_eta[i,2]/(sigma2_eta[i,2] + 1))
  tmp <- pnorm(z)
  #tmp <- pnorm( (((X.obs - mu_X[i,2])/sqrt(sigma2_eta[i,2] + sigma2_X[i,2])) - (apply(theta[,,2], 2, quantile, probs=0.5)[i])) / sqrt(sigma2_eta[i,2] / (sigma2_eta[i,2] + sigma2_X[i,2])) )
  p.obs[, i, 2] <- as.vector(tmp)
  
  
  # tmp <- pnorm((X.pred - (apply(theta[,,1], 2, quantile, probs=0.5)[i]))/ ( sqrt(apply(sigma2_eta, 2, quantile, probs=0.5)[1]) + sqrt(apply(sigma2_THETA[,,1], 2, quantile, probs=0.5)[i])))
  # p.pred[, i, 1] <- as.vector(tmp)
  # 
  # tmp <- pnorm((X.obsc - (apply(theta[,,1], 2, quantile, probs=0.5)[i]))/ ( sqrt(apply(sigma2_eta, 2, quantile, probs=0.5)[1]) + sqrt(apply(sigma2_THETA[,,1], 2, quantile, probs=0.5)[i])))
  # p.obs[, i, 1] <- as.vector(tmp)
  # 
  # 
  # tmp <- pnorm((X.pred - (apply(theta[,,2], 2, quantile, probs=0.5)[i]))/ ( sqrt(apply(sigma2_eta, 2, quantile, probs=0.5)[2]) + sqrt(apply(sigma2_THETA[,,2], 2, quantile, probs=0.5)[i])))
  # p.pred[, i, 2] <- as.vector(tmp)
  # 
  # tmp <- pnorm((X.obsc - (apply(theta[,,2], 2, quantile, probs=0.5)[i]))/ ( sqrt(apply(sigma2_eta, 2, quantile, probs=0.5)[2]) + sqrt(apply(sigma2_THETA[,,2], 2, quantile, probs=0.5)[i])))
  # p.obs[, i, 2] <- as.vector(tmp)
}





hist((dataToJags$X[dataToJags$sex==1] - mean(dataToJags$X))/sd(dataToJags$X)
     , xlim=range(X.pred)
     ,border=t_col("lightgrey",20), lty="blank"
     ,xlab="",ylab="", xaxt='n', yaxt='n',main=""
     , freq = FALSE,col = t_col("lightgrey",20)) # Save 2nd histogram data
#plot(hgB, col = c2)

# prepare graphics to add second plot
par(new = TRUE)
#par(mfrow=c(1,2))
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



hist((new.df$X[sex==2] - mean(new.df$X))/sd(new.df$X)
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




plot(NULL, xlim=c(0,4),ylim=c(0,1),xlab="Genotype", ylab="VGLL3 variances",xaxt="n")
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




plot(NULL, xlim=c(0,4),ylim=c(0,4),xlab="Genotype", ylab="Residual genetic variances",xaxt="n")
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

ratio <- samples$BUGSoutput$sims.list$ratio
plot(NULL, xlim=c(0,3),ylim=c(0,1),xlab="", ylab="Proportion residual variance",xaxt="n",main="Environment")
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

plot(NULL, xlim=c(0,3),ylim=c(0,1),xlab="", ylab="Proportion residual variance",xaxt="n",main="Genetic")
axis(1, at=c(1,2),labels=c("Male", "Female"),las=2)
#axis(1, at=1:4,labels=c("ratio_eta", "ratio_eta", "ratio_res", "ratio_res"),las=2)
#legend("topleft", c("Male", "Female"),pch=c(24,21),col=c(1,1),bty="n")
for (j in 3:4){
  segments(j-2, apply(ratio, 2, quantile, probs=0.025)[j],j-2, apply(ratio, 2, quantile, probs=0.975)[j],col="black")
  segments(j-2, apply(ratio, 2, quantile, probs=0.25)[j],j-2, apply(ratio, 2, quantile, probs=0.75)[j],col="black",lwd=2)
  points(j-2,apply(ratio, 2, quantile, probs=0.5)[j], pch=21, col="black",bg="white")
  text(j-2,apply(ratio, 2, quantile, probs=0.99)[j]+1,female,col="black")
}
  

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
par(mfrow=c(2,2))

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
a <- samples$BUGSoutput$sims.list$a
plot(NULL, xlim=c(0,3),ylim=range(a),xlab="Sex", ylab="Genotypic value",xaxt="n")
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
k <- samples$BUGSoutput$sims.list$k
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



# Contribution total variance = Percentages of the total phenotypic variance
#h2 <- samples$BUGSoutput$sims.list$h2
h2 <- mcmc[,paste0("h2[",1:2,"]")]

plot(NULL, xlim=c(0,3),ylim=c(0,1),xlab="", ylab="Total variance contribution",xaxt="n")
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
eta <- samples$BUGSoutput$sims.list$eta
#X.scaled <- (X-mean(X))/sd(X)
X.scaled=NULL
for (i in 1: length(X)){
  X.scaled[i] <- (X[i]-mu_X[dataToJags$g[i],dataToJags$sex[i]])/sigma_X[dataToJags$g[i],dataToJags$sex[i]]
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
etas <- mcmc[,paste0("eta[",1:dataToJags$N,"]")]
thetas <- mcmc[,paste0("theta[",1:dataToJags$N,"]")]

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
plot(NULL,xlim=c(-5,5),ylim=c(-5,5),xlab="Proximate cue",ylab="Threshold", main=paste0("Male (iteration: ", smp,")"))
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
#     tmp2= etas.median[dataToJags$g==gene & dataToJags$sex==sexe]
#     tmp3= thetas.median[dataToJags$g==gene & dataToJags$sex==sexe]
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

thetas <- samples$BUGSoutput$sims.list$theta
etas <- samples$BUGSoutput$sims.list$eta

#thetas_male <- thetas[,sex==1]
#thetas_female <- thetas[,sex==2]


#thetas_medians=array(,dim=c(nrow(thetas), 33,2))
thetas_medians=array(,dim=c(5, 32,2));rownames(thetas_medians)<-c("2.5%","25%","50%","75%","97.5%")
etas_medians=array(,dim=c(5, 32,2));rownames(etas_medians)<-c("2.5%","25%","50%","75%","97.5%")
etas_medians_all=array(,dim=c(5, 32));rownames(etas_medians_all)<-c("2.5%","25%","50%","75%","97.5%")

  for (y in 1987:2016){
    #thetas_medians_all[,y-1985,s] <- quantile(apply(thetas[,sex==s & year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
    etas_medians_all[,y-1985] <- quantile(apply(etas[,year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
    for (s in 1:2){
    thetas_medians[,y-1985,s] <- quantile(apply(thetas[,sex==s & year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
    etas_medians[,y-1985,s] <- quantile(apply(etas[,sex==s & year==y],1,median),probs=c(0.025,0.25, 0.5,0.75, 0.975))
  }              
}



thetas_50 <- samples$BUGSoutput$median$theta
etas_50 <- samples$BUGSoutput$median$eta
thetas_medians_ind=etas_medians_ind=array(,dim=c(dataToJags$N,32,2))
for (s in 1:2){
  for (y in 1987:2016){
    tmp <- thetas_50[sex==s & year==y]
    thetas_medians_ind[1:length(tmp),y-1985,s] <- tmp
    tmp=NULL
    tmp <- etas_50[sex==s & year==y]
    etas_medians_ind[1:length(tmp),y-1985,s] <- tmp
  }
}


#par(mfcol=c(1,2))

# #Xc <- (X-mean(X))/sd(X)
# X.scaled=NULL
# for (i in 1: length(X)){
#   X.scaled[i] <- (X[i]-mu_X[dataToJags$g[i],dataToJags$sex[i]])/sigma_X[dataToJags$g[i],dataToJags$sex[i]]
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

range_years <- 1986:2017
plot(NULL, xlim=c(1986,2017),ylim=c(-3,3), ylab="Proximate cue",xlab="")
segments(range_years-.1,etas_medians_all["2.5%",], range_years-.1,etas_medians_all["97.5%",] ,col=1)
segments(range_years-.1,etas_medians_all["25%",], range_years-.1,etas_medians_all["75%",] ,col=1,lwd=2)
points(range_years-.1,etas_medians_all["50%",],col=1,pch=16)
# for (y in 1:32){
#   for (i in 1:nrow(etas_medians_ind[,,1])){
#  #   points(1985+y, etas_medians_ind[i,y,1],col=1,pch=3)
#   }
#   }

y <- etas_medians_all["50%",]
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
plot(NULL, xlim=c(1986,2017),ylim=c(-3,3), ylab="Tresholds",xlab="")

# for (y in 1986:2017){
# points(y-.1,thetas_50[year==y & ],col=1,pch=16)
# }

segments(range_years-.1,thetas_medians["2.5%",,1], range_years-.1,thetas_medians["97.5%",,1] ,col=1)
segments(range_years-.1,thetas_medians["25%",,1], range_years-.1,thetas_medians["75%",,1] ,col=1,lwd=2)
points(range_years-.1,thetas_medians["50%",,1],col=1,pch=16)
# for (y in 1:32){
#   for (i in 1:nrow(thetas_medians_ind[,,1])){
# #    points(1985+y, thetas_medians_ind[i,y,1],col=1,pch=3)
#   }
# }

y <- thetas_medians["50%",,1]
x <- range_years
mod1=loess(y~x,span=0.5)
xfit=seq(from=min(x),to=max(x),length.out=30)
yfit1=predict(mod1,newdata=xfit)
points(xfit,yfit1,type="l",lwd=2,col=1)


segments(range_years+.1,thetas_medians["2.5%",,2], range_years+.1,thetas_medians["97.5%",,2] ,col=2)
segments(range_years+.1,thetas_medians["25%",,2], range_years+.1,thetas_medians["75%",,2] ,col=2,lwd=2)
points(range_years+.1,thetas_medians["50%",,2],col=2,pch=16)
for (y in 1:32){
  for (i in 1:nrow(thetas_medians_ind[,,2])){
 #   points(1985+y, thetas_medians_ind[i,y,2],col=2,pch=3)
  }
}

y <- thetas_medians["50%",,2]
x <- range_years
mod1=loess(y~x,span=0.5)
xfit=seq(from=min(x),to=max(x),length.out=30)
yfit1=predict(mod1,newdata=xfit)
points(xfit,yfit1,type="l",lwd=2,col=2)


#linearMod <- lm(y ~ x)  # build linear regression model on full data
#print(linearMod)


plot(thetas_medians["50%",,1],thetas_medians["50%",,2],pch="", xlab="MALE",ylab="FEMALE", main="Tresholds (medians)")
text(thetas_medians["50%",,1],thetas_medians["50%",,2], 1986:2017)
cor.test( thetas_medians["50%",,1],thetas_medians["50%",,2])

#### ALLELE FREQUENCIES
#par(mfcol=c(1,1))
tmp <- aggregate(X~g+sex+t, data=new.df, FUN=length)
n <- xtabs( X ~ g+sex+t, tmp) # convert dataframe to matrix
## Frequence of allele E
p=q=array(,dim=c(2,32))
for (s in 1:2){    # sex
  for (t in 2:31){
    p[s,t] <- sum((2*n[1,s,t])+(1*n[2,s,t])+(0*n[3,s,t]))/(2*sum(n[1:3,s,t])) # population frequence allelic for E
    q[s,t] <- 1 - p[s,t]
  }}
colnames(p)<-c(range_years);rownames(p)<-c("Male","Female")

#par(mfcol=c(1,1))
plot(NULL,xlim=c(1986,2017),ylim=c(0.5,1),xlab="",ylab="Allele E frequencies")
for (s in 1:2){
  y <- p[s,]
  x <- range_years
  mod1=loess(y~x,span=0.9)
  xfit=seq(from=min(x),to=max(x),length.out=30)
  yfit1=predict(mod1,newdata=xfit)
  points(xfit,yfit1,type="l",lwd=2,col=col.sex[s])
  points(range_years,p[s,],col=col.sex[s],pch=16)
}
legend("topleft", legend=c("Male", "Female"),fill=col.sex, bty="n",border = col.sex, title ="Allele E frequencies")




par(mfcol=c(2,2))
plot(etas_medians["50%",2:31,1],thetas_medians["50%",2:31,1], xlim=c(-1,1), ylim=c(-2,-0.5), pch="",main="Male", xlab="Proximate cues",ylab="Thresholds")
points(etas_medians["50%",2:31,1],thetas_medians["50%",2:31,1], pch="")
text(etas_medians["50%",2:31,1],thetas_medians["50%",2:31,1],1986:2016,cex=0.75)

plot(etas_medians["50%",2:31,2],thetas_medians["50%",2:31,2], xlim=c(-.5,1), ylim=c(0.25,1.25), pch="", main="Female", xlab="Proximate cues",ylab="Thresholds")
points(etas_medians["50%",2:31,2],thetas_medians["50%",2:31,2], pch="")
text(etas_medians["50%",2:31,2],thetas_medians["50%",2:31,2],1986:2016,cex=0.75)


#par(mfcol=c(1,2))
plot(p[1,2:31],thetas_medians["50%",2:31,1], xlim=c(0.5,1), ylim=c(-2,-0.5), pch="",main="Male", xlab="Freq allele E",ylab="Thresholds")
points(p[1,2:31],thetas_medians["50%",2:31,1], pch="")
text(p[1,2:31],thetas_medians["50%",2:31,1],1986:2016,cex=0.75)

plot(p[2,2:31],thetas_medians["50%",2:31,2], xlim=c(0.5,1), ylim=c(0.25,1.25), pch="", main="Female", xlab="Freq allele E",ylab="Thresholds")
points(p[2,2:31],thetas_medians["50%",2:31,2], pch="")
text(p[2,2:31],thetas_medians["50%",2:31,2],1986:2016,cex=0.75)








par(mfcol=c(1,1))

thetas <- samples$BUGSoutput$sims.list$theta
etas <- samples$BUGSoutput$sims.list$eta


thetas.stat <- apply(thetas,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
etas.stat <- apply(etas,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

smp=sample(1:dataToJags$N,50,replace=FALSE)

plot(NULL, xlim=c(0,50), ylim=c(-10,10),xlab="",ylab="Individuals")
for (i in 1:50){
  segments(i-.1,thetas.stat["2.5%",smp[i]], i-.1,thetas.stat["97.5%",smp[i]] ,col=1)
  segments(i-.1,thetas.stat["25%",smp[i]], i-.1,thetas.stat["75%",smp[i]] ,col=1,lwd=2)
  points(i-0.1, thetas.stat[3,smp[i]], pch=21,bg="white",col=1)
  
  segments(i+.1,etas.stat["2.5%",smp[i]], i+.1,etas.stat["97.5%",smp[i]] ,col=2)
  segments(i+.1,etas.stat["25%",smp[i]], i+.1,etas.stat["75%",smp[i]] ,col=2,lwd=2)
  points(i+.1, etas.stat[3,smp[i]], pch=21,bg="white",col=2)
  
  #points(i+.1, dataToJags$X.scaled[smp[i]], pch="X",bg="white",col=2)
  points(i+.1, X.scaled[smp[i]], pch="X",bg="white",col=2)
  text(i,-10,labels=dataToJags$Y[smp[i]],cex=0.75)  
}
legend("topleft",legend=c("Threshold","Prox. Cue","Growth scaled"),pch=c(1,1,"X"),col=c(1,2,2) ,bty="n")


dev.off()


