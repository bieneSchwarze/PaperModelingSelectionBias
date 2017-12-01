########################################################################
########################################################################
##
## SC4 - math competence
##
## Analyses : Linear Model with Correction for Selection Bias
##
## Method: MI (regarding the multilevel structure when imputing values)
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
library(mice)                                                                
library(miceadds)
# install: devtools::install_github("bbolker/mice")
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
# LITTLE's MCAR test
# --------------------------------------------------------------------
datM <- dat[,c("ID_t", "ID_i", "sc", "dgcf", "math1", "math2")]    
littleTest <- LittleMCAR(datM[,-which(names(datM) %in% c("ID_t", "ID_i"))])
littleTest$p.value
# small p-value: much evidence against MCAR (H_0)

# --------------------------------------------------------------------
# MI Model 
# --------------------------------------------------------------------

# Reshape data to long format
dat$sc <- dat$sc-mean(dat$sc)
dat$sc <- dat$sc/sqrt(var(dat$sc))
dat$dgcf <- dat$dgcf-mean(dat$dgcf)
dat$dgcf <- dat$dgcf/sqrt(var(dat$dgcf))
datLongM <- reshape(dat, idvar="ID_t",
                    varying=list(math=c("math1", "math2")),
                    v.names=c("math"),direction="long", sep="")
datLongM <- datLongM[order(datLongM[,1]),]
datLongM$time <- ifelse(datLongM$time == 2, 1, 0)
datLongM <- as.numeric.matrix(datLongM)

ini <- mice(datLongM,maxit=0)
pred <- ini$pred
pred["math","ID_i"] <- 0
pred["math", "ID_t"] <- -2
meth <- ini$meth
meth["math"] <- "2l.pan" 
miR <- 20
imp <- mice(datLongM, m=miR, pred=pred, meth=meth, seed=1012)  

# Pool results
estMI <- function(){
  cat("It: ",cc, "\n\n")
  datMI.long <- complete(imp, action=cc, include=TRUE)
  cc <<- cc + 1
  datMI.long <- datMI.long[order(datMI.long[,1]),]  
  
  sc <- datMI.long[!duplicated(datMI.long[,1]),"sc"]  
  scM <- sc - mean(sc) 
  scM <- scM/sqrt(var(scM))
  scI <- cbind(datMI.long[!duplicated(datMI.long[,1]),1], scM)
  colnames(scI) <- c("ID_t", "scM")
  datMI.long <- merge(datMI.long, scI, by="ID_t")
  datMI.long <- datMI.long[order(datMI.long[,1]),]   
  
  dgcf <- datMI.long[!duplicated(datMI.long[,1]),"dgcf"]  
  dgcfM <- dgcf - mean(dgcf) 
  dgcfM <- dgcfM/sqrt(var(dgcfM))
  dgcfI <- cbind(datMI.long[!duplicated(datMI.long[,1]),1], dgcfM)
  colnames(dgcfI) <- c("ID_t", "dgcfM")
  datMI.long <- merge(datMI.long, dgcfI, by="ID_t")
  datMI.long <- datMI.long[order(datMI.long[,1]),]      

  res <- lmer(math ~ scM + dgcfM + time*scM + time*dgcfM + (1|ID_i) + (1|ID_t), 
          data=datMI.long)
  return(res)
}
cc <- 1
fitP <- with(imp, estMI())
pooled.mi <- pool(fitP)
fitRes <- summary(pooled.mi)

# Pick random effects
rf <- matrix(NA, nrow=miR, ncol=3)
for(i in 1:miR){
  smod <- summary(fitP[[4]][[i]])
  rf[i,] <- c(unlist(smod$varcor),smod$sigma^2)
}
ranEffMI <- apply(rf, 2, median)
 
# ------------------------------------------------------------------------
# DERIVE CONFIDENCE INTERVALS FOR RANDOM EFFECTS VIA CLUSTER BOOTSTRAPPING
# ------------------------------------------------------------------------

nr <- 2
nCL <- length(unique(dat[,1]))
allCL <- unique(dat[,1])
seedList <- sample(x=10000, size=nr)

doTheVarMI <- function(it){

  datV <- dat[dat[,1] %in% sample(allCL, nCL-1, replace=TRUE),]      
  
  datLongV <- reshape(datV, idvar="ID_t",
                    varying=list(math=c("math1", "math2")),
                    v.names=c("math"),direction="long", sep="")
  datLongV <- datLongV[order(datLongV[,1]),]
  datLongV$time <- ifelse(datLongV$time == 2, 1, 0)
  datLongV <- as.numeric.matrix(datLongV)   

  ini <- mice(datLongV,maxit=0)
  pred <- ini$pred
  pred["math","ID_i"] <- 0
  pred["math", "ID_t"] <- -2
  meth <- ini$meth
  meth["math"] <- "2l.pan"
  miR <- 20
  imp <- mice(data=datLongV,
              m=miR,
              imputationMethod = meth,
              predictorMatrix=pred,
              seed=1012)     
  
  estMI <- function(){

    datMI.long <- complete(imp, action=cc, include=TRUE)
    cc <<- cc + 1
    datMI.long <- datMI.long[order(datMI.long[,1]),]  
    
    sc <- datMI.long[!duplicated(datMI.long[,1]),"sc"]  
    scM <- sc - mean(sc) 
    scM <- scM/sqrt(var(scM))
    scI <- cbind(datMI.long[!duplicated(datMI.long[,1]),1], scM)
    colnames(scI) <- c("ID_t", "scM")
    datMI.long <- merge(datMI.long, scI, by="ID_t")
    datMI.long <- datMI.long[order(datMI.long[,1]),]   
    
    dgcf <- datMI.long[!duplicated(datMI.long[,1]),"dgcf"]  
    dgcfM <- dgcf - mean(dgcf) 
    dgcfM <- dgcfM/sqrt(var(dgcfM))
    dgcfI <- cbind(datMI.long[!duplicated(datMI.long[,1]),1], dgcfM)
    colnames(dgcfI) <- c("ID_t", "dgcfM")
    datMI.long <- merge(datMI.long, dgcfI, by="ID_t")
    datMI.long <- datMI.long[order(datMI.long[,1]),]         

    res <- lmer(math ~ scM + dgcfM + time*scM + time*dgcfM + (1|ID_i) + (1|ID_t), 
            data=datMI.long)    
      
    return(res)
  }
  cc <- 1
  fitP <- with(imp, estMI())
  pooled.mi <- pool(fitP)
  fitRes <- summary(pooled.mi)

  rf <- matrix(NA, nrow=miR, ncol=3)
  for(i in 1:miR){
    smod <- summary(fitP[[4]][[i]])
    rf[i,] <- c(unlist(smod$varcor),smod$sigma^2)
  }
  ranEffMI <- apply(rf, 2, median)
    
  return(list(coeffSelB=fitRes[,1], ranEffMI=ranEffMI)) 
}

sfInit(parallel=T, cpus=8, type='SOCK', slaveOutfile="varJK_PARALLEL_MI_20170818.log")
sfLibrary(lme4)
sfLibrary(reshape)
sfLibrary(mice)
sfLibrary(rms)
sfExportAll()
res.v <- sfLapply(1:nr, doTheVarMI)  # bootstrap
sfStop()

contMIBeta <- matrix(NA, nrow=nr, ncol=length(names(res.v[[1]]$coeffSelB)))
colnames(contMIBeta) <- names(res.v[[1]]$coeffSelB)
contMIRand <-  matrix(NA, nrow=nr, ncol=3)

for(i in 1:nr){
   liEl <- res.v[[i]]
   r <- liEl$coeffSelB
   contMIBeta[i,which(colnames(contMIBeta)%in% names(r))] <- r
   contMIRand[i,] <- c(liEl$ranEffMI)
}

# Remove lines with NAs   (occurs if school types are not part of the bootstrap sample, e.g. RS)
contMIBeta_n <- contMIBeta
rem <- which(is.na(contMIBeta_n), arr.ind=TRUE)[,1]
if(length(rem)>0){contMIBeta <- contMIBeta_n[-rem,]; contMIRand <- contMIRand[-rem,]}

colnames(contMIRand) <- c("sig.id", "sig.inst", "residual") 

# There are several variables for which bootstrap samples are not centered and symmetric. 
estimMI <- c(fitRes[,1], ranEffMI)
bootM_MI <- cbind(contMIBeta, contMIRand) 
tabBOOT_MI <- bootBCa(estimMI, bootM_MI, type="basic", n=nrow(contMIBeta), seed=123, conf.int = 0.95) 

CI.clusterMI <- cbind(tabBOOT_MI[1,], estimMI, tabBOOT_MI[2,])
colnames(CI.clusterMI) <- c("lower CI", "estim", "upper CI")
CI.clusterMI[1:6,1] <-  fitRes[1:6,6] 
CI.clusterMI[1:6,3] <-  fitRes[1:6,7] 
