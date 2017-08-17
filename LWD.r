########################################################################
########################################################################
##
## SC4 - math competence
##
## Analyses : Growth Curve Model (Linear panel model) 
##
## Method: Listwise deletion 
##
## @project Modeling longitudinal dropout
## @author Sabine Zinn
##
## 16.08.2017
##
########################################################################
########################################################################

rm(list=ls())

# load relevant packages
library(lme4)
library(reshape)
library(snowfall)
library(rms)

# --------------------------------------------------------------------
# LOAD DATA
# --------------------------------------------------------------------

setwd("Z:\\Projektgruppen_(08)\\StatMethoden_(p000033)\\SelectivityLinPanelModel\\ModelCodeVersion1\\syntax\\SupplementalMaterialPaper\\") 
dat <- read.table("dataSC4MatheWIDE.csv", sep=",", na.strings = ".")
colnames(dat) <- c("ID_t","ID_i","dgcf","sc","sex","HS","RS","REST","migration","amode","age","math1","math2","DROP")
                                 
# aux functions
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
as.numeric.matrix <- function(mat){
 for(i in 1:ncol(mat)){
    if(is.factor(mat[,i]))
      mat[,i] <- as.numeric.factor(mat[,i])
 }
 return(as.matrix(mat))
}

# --------------------------------------------------------------------
# RESHAPE RELEVANT VARIABLES TO LONG FORMAT
# --------------------------------------------------------------------

datM <- dat[,c("ID_t", "ID_i", "sc", "dgcf", "math1", "math2")]
datRed <-  datM[!is.na(datM$math2),]  
datRed$sc <- datRed$sc-mean(datRed$sc)
datRed$sc <- datRed$sc/sqrt(var(datRed$sc))
datRed$dgcf <- datRed$dgcf-mean(datRed$dgcf)
datRed$dgcf <- datRed$dgcf/sqrt(var(datRed$dgcf))  
datLongM <- reshape(datRed, idvar="ID_t",
                    varying=list(math=c("math1", "math2")),
                    v.names=c("math"),direction="long", sep="")
datLongM <- datLongM[order(datLongM[,1]),]
datLongM$time <- ifelse(datLongM$time == 2, 1, 0)
datLongM <- as.numeric.matrix(datLongM)


# --------------------------------------------------------------------
# ESTIMATE MODEL & DERIVE CONFIDENCE INTERVALS FOR EFFECT ESTIMATES
# --------------------------------------------------------------------

anaModNull <- lmer(math ~ sc + dgcf + time*sc + time*dgcf + (1|ID_i) + (1|ID_t), data=as.data.frame(datLongM))
sumAnaNull <- summary(anaModNull)
coefAnaMod <- coef(summary(sumAnaNull))[,1]
stdCoefAnaMod <- coef(summary(sumAnaNull))[,2]
alpha <- 0.05
ciL <- coefAnaMod - 1.96*qt(1-alpha/2,nrow(datLongM))*stdCoefAnaMod
ciU <- coefAnaMod + 1.96*qt(1-alpha/2,nrow(datLongM))*stdCoefAnaMod
CI.clusterNull <- cbind(ciL, coefAnaMod, ciU)
ranEffNull <- unlist(VarCorr(anaModNull))
resNull <- sigma(anaModNull)^2

# ------------------------------------------------------------------------
# DERIVE CONFIDENCE INTERVALS FOR RANDOM EFFECTS VIA CLUSTER BOOTSTRAPPING
# ------------------------------------------------------------------------

nr <- 200
nCL <- length(unique(datRed[,1]))
allCL <- unique(datRed[,1])
seedList <- sample(x=10000, size=nr)

doTheVar <- function(it){

  datV <- datRed[datRed[,1] %in% sample(allCL, nCL-1, replace=TRUE),]      
  datV$sc <- datV$sc-mean(datV$sc)
  datV$sc <- datV$sc/sqrt(var(datV$sc))
  datV$dgcf <- datV$dgcf-mean(datV$dgcf)
  datV$dgcf <- datV$dgcf/sqrt(var(datV$dgcf))  
  
  datLongV <- reshape(datV, idvar="ID_t",
                    varying=list(math=c("math1", "math2")),
                    v.names=c("math"),direction="long", sep="")
  datLongV <- datLongV[order(datLongV[,1]),]
  datLongV$time <- ifelse(datLongV$time == 2, 1, 0)
  datLongV <- as.numeric.matrix(datLongV)   

  anaModNull <- lmer(math ~ sc + dgcf + time*sc + time*dgcf + (1|ID_i) + (1|ID_t), data=as.data.frame(datLongV))
  sumAnaNull <- summary(anaModNull)
  coefAnaMod <- coef(summary(sumAnaNull))[,1]
  stdCoefAnaMod <- coef(summary(sumAnaNull))[,2]
  ciL <- coefAnaMod - 1.96*qt(1-alpha/2,nrow(datLongV))*stdCoefAnaMod
  ciU <- coefAnaMod + 1.96*qt(1-alpha/2,nrow(datLongV))*stdCoefAnaMod
  CI.clusterNull <- cbind(ciL, coefAnaMod, ciU)
  ranEffNull <- unlist(VarCorr(anaModNull))
  resNull <- sigma(anaModNull)^2
    
  return(list(coeffSelB=coefAnaMod, ranEff=c(ranEffNull, resNull))) 
}

sfInit(parallel=T, cpus=6, type='SOCK', slaveOutfile="varJK_PARALLEL_LWD_20170818.log")
sfLibrary(lme4)
sfLibrary(nnet)
sfLibrary(reshape)
sfExportAll()
res.v <- sfLapply(1:nr, doTheVar)  # bootstrap
sfStop()

contBeta <- matrix(NA, nrow=nr, ncol=length(names(res.v[[1]]$coeffSelB)))
colnames(contBeta) <- names(res.v[[1]]$coeffSelB)
contRand <-  matrix(NA, nrow=nr, ncol=3)

for(i in 1:nr){
   liEl <- res.v[[i]]
   r <- liEl$coeffSelB
   contBeta[i,which(colnames(contBeta)%in% names(r))] <- r
   contRand[i,] <- c(liEl$ranEff)
}

# Remove lines with NAs   (occurs if school types are not part of the bootstrap sample, e.g. RS)
contBeta_n <- contBeta
rem <- which(is.na(contBeta_n), arr.ind=TRUE)[,1]
if(length(rem)>0){contBeta <- contBeta_n[-rem,]; contMIRand <- contMIRand[-rem,]}

colnames(contRand) <- c("sig.id", "sig.inst", "residual") 

# There are several variables for which bootstrap samples are not centered and symmetric. 
estim <- c(coefAnaMod, c(ranEffNull, resNull))
bootM <- cbind(contBeta, contRand) 
tabBOOT <- bootBCa(estim, bootM, type="basic", n=nrow(contBeta), seed=123, conf.int = 0.95) 

CI.cluster <- cbind(tabBOOT[1,], estim, tabBOOT[2,])
colnames(CI.cluster) <- c("lower CI", "estim", "upper CI")
CI.cluster[1:6,] <-  CI.clusterNull 
