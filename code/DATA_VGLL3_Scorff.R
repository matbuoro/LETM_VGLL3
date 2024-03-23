
#________________________DATA________________________#
#### DATA -------
load("data/data_vgll3+year.rdata")
df <- df[order(df$t),]
#View(df)

new.df <- na.omit(df)

# remove year 1985 (only 1 female)
new.df <- new.df[-which(new.df$t==1985),]


#attach(new.df)

# g=as.numeric(new.df$g)
# sex=as.numeric(new.df$sex)
# year = new.df$t




# Number of fish by genotype and sex:
new.df$period <- ifelse(new.df$t<2005,1,2) # after 2004
tmp <- aggregate(X~g+sex+period, data=new.df, FUN=length)
n <- xtabs( X ~ g+sex+period, tmp) # convert dataframe to matrix

## Frequence of allele E
p=q=array(,dim=c(2,2))
for (s in 1:2){    # sex
  for (t in 1:2){
    p[s,t] <- sum((2*n[1,s,t])+(1*n[2,s,t])+(0*n[3,s,t]))/(2*sum(n[1:3,s,t])) # population frequence allelic for E
    q[s,t] <- 1 - p[s,t]
  }}
colnames(p)<-c("Before2006","After2005");rownames(p)<-c("Male","Female")


#sum(nByYear[,1:21])
#sum(nByYear[,22:33])
#X <- new.df$X

# Number of fish by genotype and sex:
tmp <- aggregate(X~g+sex, data=new.df, FUN=length)
n <- xtabs( X ~ g+sex, tmp) # convert dataframe to matrix


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
  , Y=new.df$Y # maturation decision: 1 for 1SW, 0 otherwise
  , X=new.df$X # Observable cue, i.e growth centre/reduit
  , X.scaled=X.scaled # Observable cue, i.e growth centre/reduit
  , mean.X=mean.X, sd.X=sd.X
  , sex=as.numeric(new.df$sex) # sex
  , g=as.numeric(new.df$g) # genotypes: 1 for EE, 2 for EL and 3 for LL
  , n=n # number of fish per genotype (row) and sex (col)
  , year = new.df$t# - (min(new.df$t)-1)
  , period = new.df$period # after 2004
  #,p=p,q=1-p
)
#attach(dataToJags)
## 
# new.df$X.c <- new.df$X - mean(new.df$X)
# period = ifelse(new.df$t<=2005,1,2)
# boxplot(X.c~period+g, data=new.df); abline(h=0,lty=2,lwd=2)