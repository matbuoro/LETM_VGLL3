
#________________________DATA________________________#
#### DATA -------
#load("data/data_vgll3+year.rdata")
#df <- df[order(df$t),]
#View(df)
# Save data frame to a text file with comments
#write.csv(df, "data/Data_vgll3_Scorff.csv",row.names = FALSE)
#df <- read.csv("data/Data_vgll3_Scorff.csv")
load("data/data_vgll3.rdata")
new.df <- na.omit(df)

# remove year 1985 (only 1 female) and last year (no MSW)
new.df <- subset(new.df, !(t %in% c(1985, 1986, 1992, 2017)))


#attach(new.df)

# g=as.numeric(new.df$g)
# sex=as.numeric(new.df$sex)
new.df$year = new.df$t
new.df$Y = new.df$Y + 1 # 1 for 1SW, 2 for MSW



# Number of fish by genotype and sex:
new.df$period <- ifelse(new.df$t<2005,1,2) # after 2004
#tmp <- aggregate(X~g+sex+period, data=new.df, FUN=length)
##n <- xtabs( X ~ g+sex+period, tmp) # convert dataframe to matrix

tmp <- aggregate(X~g+sex+Y, data=new.df, FUN=length)
counts <- xtabs( X ~ g+sex+Y, tmp) # convert dataframe to matrix

## Frequence of allele E
#p=q=array(,dim=c(2,2))
#for (s in 1:2){    # sex
#  for (t in 1:2){
#    p[s,t] <- sum((2*n[1,s,t])+(1*n[2,s,t])+(0*n[3,s,t]))/(2*sum(n[1:3,s,t])) # population frequence allelic for E
#    q[s,t] <- 1 - p[s,t]
#  }}
#colnames(p)<-c("Before2006","After2005");rownames(p)<-c("Male","Female")


# Convert xtabs to array
counts <- aperm(array(counts, dim = c(3, 2, 2)), c(1, 2, 3))

freqs <- array(,dim=c(3,2,2))
for (i in 1:2) {
  for (j in 1:2) {
    for (k in 1:3){
      freqs[k,i,j] <- counts[k,i,j]/sum(counts[1:3,i,j])
    }
  }
}

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


mat <- array(,dim=c(3,2,2))
mat[,1,1]<- freqs[1:3,1,1]* p_1SW * p_male_1SW
mat[,2,1] <- freqs[1:3,2,1]* p_1SW * p_female_1SW
mat[,1,2]<- freqs[1:3,1,2]* (1-p_1SW) * p_male_MSW
mat[,2,2] <- freqs[1:3,2,2]* (1-p_1SW) * p_female_MSW

sum(mat) # should be 1!



#sum(nByYear[,1:21])
#sum(nByYear[,22:33])
#X <- new.df$X

# Number of fish by genotype and sex:
#tmp <- aggregate(X~g+sex, data=new.df, FUN=length)
#n <- xtabs( X ~ g+sex, tmp) # convert dataframe to matrix


tmp <- aggregate(X~g+sex, data=new.df, FUN=mean)
mean.X <- xtabs( X ~ g+sex, tmp) # convert dataframe to matrix
tmp <- aggregate(X~g+sex, data=new.df, FUN=sd)
sd.X <- xtabs( X ~ g+sex, tmp) # convert dataframe to matrix
# gr=NULL
# g=as.numeric(new.df$g)
# for (i in 1:nrow(new.df)) {
#   if(new.df$sex[i]==1) gr[i] <- g[i]
#   if(new.df$sex[i]==2) gr[i] <- g[i] + 3
# }
# pairwise.t.test(X, gr, pool.sd = FALSE, paired = FALSE,alternative = c("two.sided"))

X.scaled=NULL
for (i in 1:nrow(new.df)){
  X.scaled[i] <- (new.df$X[i]-mean.X[new.df$g[i],new.df$sex[i]])/sd.X[new.df$g[i],new.df$sex[i]] # scaled growth
}

# tmp <- aggregate(X.scaled~g+sex, FUN=mean)
# mean.X <- xtabs( X.scaled ~g+sex, tmp) # convert dataframe to matrix
# tmp <- aggregate(X.scaled~g+sex, FUN=sd)
# sd.X <- xtabs( X.scaled~g+sex, tmp) # convert dataframe to matrix

#X.scaled <- (X-mean(X))/sd(X) #scale(new.df$X, center = TRUE, scale = TRUE))

## Data used for inference in jags
dataToJags <- list(
  N=nrow(new.df) # total number of fish
  , Y=new.df$Y-1 # maturation decision: 1 for 1SW, 0 otherwise
  , X=new.df$X # Observable cue, i.e growth centre/reduit
  , X.scaled=X.scaled # Observable cue, i.e growth centre/reduit
  , mean.X=mean.X, sd.X=sd.X
  , sex=as.numeric(new.df$sex) # sex
  , g=as.numeric(new.df$g) # genotypes: 1 for EE, 2 for EL and 3 for LL
  #, n=n # number of fish per genotype (row) and sex (col)
  , year = new.df$t# - (min(new.df$t)-1)
  , period = new.df$period # after 2004
  ,freq = mat[,,1]+mat[,,2]
  #,p=p,q=1-p
)
#attach(dataToJags)
## 
# new.df$X.c <- new.df$X - mean(new.df$X)
# period = ifelse(new.df$t<=2005,1,2)
# boxplot(X.c~period+g, data=new.df); abline(h=0,lty=2,lwd=2)
