########################################################################
########################################################################
##
## NEPS Starting Cohort 4 - math competence
##
## Prepare data for further analyses
##
## @project Modeling longitudinal dropout
## @author Timo Gnambs & Sabine Zinn
##
## 16.08.2017
##
########################################################################
########################################################################

rm(list = ls())

# load relevant packages
library(haven)
library(doBy)
library(MBESS)
library(lavaan)
library(mice)

# --------------------------------------------------------------------
# LOAD DATA
# --------------------------------------------------------------------

# set path
path <- "Z:\\Projektgruppen_(08)\\StatMethoden_(p000033)\\SelectivityLinPanelModel\\rawdata"

# set paths to SUFs
path.SC4 <- paste(path, "SC4 9.0.0", sep = "\\")

# load SUFs
wd <- getwd()
setwd(path.SC4)
dat1  <- read_spss("SC4_pTarget_D_9-0-0.sav")
dat2  <- read_spss("SC4_xTargetCompetencies_D_9-0-0.sav")
dat3  <- read_spss("SC4_pParent_D_9-0-0.sav")
dat4  <- read_spss("SC4_CohortProfile_D_9-0-0.sav")
dat5  <- read_spss("SC4_Weights_D_9-0-0.sav")  
dat6  <- read.table("SC4_Studynumber_w7.txt", sep="\t", header=TRUE)
setwd(wd)
rm(wd, path.SC4)

# --------------------------------------------------------------------
# MERGE FILES
# --------------------------------------------------------------------

# select competences
sel <- c("ID_t", "wave_w1", "wave_w2", "wave_w7",
         # DGCF
         "dgci2101_sc4g9_c", "dgci2102_sc4g9_c", "dgci2103_sc4g9_c",
         "dgci2104_sc4g9_c", "dgci2201_sc4g9_c", "dgci2202_sc4g9_c",
         "dgci2203_sc4g9_c", "dgci2204_sc4g9_c", "dgci2301_sc4g9_c",
         "dgci2302_sc4g9_c", "dgci2303_sc4g9_c", "dgci2304_sc4g9_c",
         "dgg9_sc3b",
         # mathematical competence (WLE)
         "mag9_sc1", "mag12_sc1u")
sc4 <- dat2[, sel]
rm(sel)

# merge with self-concept
sc4 <- merge(sc4,
             dat1[dat1$wave == 1,
                  c("ID_t", "t66001a", "t66001b", "t66001c")],
             by = "ID_t", all.x = TRUE)

# merge with sociodemographics
sc4 <- merge(sc4,
             dat1[dat1$wave == 1,
                  c("ID_t", "t70004m", "t70004y", # birth year / month
                    "t700031", "t400500_g1")],    # sex, migration
             by = "ID_t", all.x = TRUE)
sc4 <- merge(sc4,
             dat1[dat1$wave == 2,
                  c("ID_t", "t70004m", "t70004y", # birth year / month
                    "t700031", "t400500_g1")],    # sex, migration
             by = "ID_t", all.x = TRUE, suffixes = c(".w1", ".w2"))
sc4 <- merge(sc4,
             dat1[dat1$wave == 3,
                  c("ID_t", "t70004m", "t70004y", # birth year / month
                    "t700031", "t400500_g1")],    # sex, migration
             by = "ID_t", all.x = TRUE)
sc4 <- merge(sc4,
             dat1[dat1$wave == 4,
                  c("ID_t", "t70004m", "t70004y", # birth year / month
                    "t700031", "t400500_g1")],    # sex, migration
             by = "ID_t", all.x = TRUE, suffixes = c(".w3", ".w4"))
sc4 <- merge(sc4,
             dat1[dat1$wave == 5,
                  c("ID_t", "t70004m", "t70004y", # birth year / month
                    "t700031", "t400500_g1")],    # sex, migration
             by = "ID_t", all.x = TRUE)
sc4 <- merge(sc4,
             dat1[dat1$wave == 6,
                  c("ID_t", "t70004m", "t70004y", # birth year / month
                    "t700031", "t400500_g1")],    # sex, migration
             by = "ID_t", all.x = TRUE, suffixes = c(".w5", ".w6"))

# merge with institution ID
sc4 <- merge(sc4, dat4[dat4$wave == 1, c("ID_t", "ID_i")],
             by = "ID_t", all.x = TRUE)

# schooltype (gen.)
sc4 <- merge(sc4, dat4[dat4$wave == 1, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE)
sc4 <- merge(sc4, dat4[dat4$wave == 2, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE, , suffixes = c(".w1", ".w2"))
sc4 <- merge(sc4, dat4[dat4$wave == 3, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE, , suffixes = c(".w2", ".w3"))
sc4 <- merge(sc4, dat4[dat4$wave == 4, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE, , suffixes = c(".w3", ".w4"))
sc4 <- merge(sc4, dat4[dat4$wave == 5, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE, , suffixes = c(".w4", ".w5"))
sc4 <- merge(sc4, dat4[dat4$wave == 6, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE, , suffixes = c(".w5", ".w6"))
sc4 <- merge(sc4, dat4[dat4$wave == 7, c("ID_t", "t723080_g1")],
             by = "ID_t", all.x = TRUE, , suffixes = c(".w6", ".w7"))
sc4 <- merge(sc4, dat5[, c("ID_t", "stratum_exp")], 
             by="ID_t", all.x=TRUE)             

# --------------------------------------------------------------------
# RECODE VARIABLES
# --------------------------------------------------------------------

# sex (0 = male, 1 = female)
sc4$sex                 <- sc4$t700031.w1
sc4$sex[is.na(sc4$sex)] <- sc4$t700031.w2[is.na(sc4$sex)]
sc4$sex[is.na(sc4$sex)] <- sc4$t700031.w3[is.na(sc4$sex)]
sc4$sex[is.na(sc4$sex)] <- sc4$t700031.w4[is.na(sc4$sex)]
sc4$sex[is.na(sc4$sex)] <- sc4$t700031.w5[is.na(sc4$sex)]
sc4$sex[is.na(sc4$sex)] <- sc4$t700031.w6[is.na(sc4$sex)]
sc4$sex <- recodeVar(sc4$sex, c(1, 2), c(0, 1))

# age (in years via century month code)
sc4$birthY <-                 sc4$t70004y.w1
sc4$birthY[is.na(sc4$birthY)] <- sc4$t70004y.w2[is.na(sc4$birthY)]
sc4$birthY[is.na(sc4$birthY)] <- sc4$t70004y.w3[is.na(sc4$birthY)]
sc4$birthY[is.na(sc4$birthY)] <- sc4$t70004y.w4[is.na(sc4$birthY)]
sc4$birthY[is.na(sc4$birthY)] <- sc4$t70004y.w5[is.na(sc4$birthY)]
sc4$birthY[is.na(sc4$birthY)] <- sc4$t70004y.w6[is.na(sc4$birthY)]
sc4$birthM <-                 sc4$t70004m.w1
sc4$birthM[is.na(sc4$birthM)] <- sc4$t70004m.w2[is.na(sc4$birthM)]
sc4$birthM[is.na(sc4$birthM)] <- sc4$t70004m.w3[is.na(sc4$birthM)]
sc4$birthM[is.na(sc4$birthM)] <- sc4$t70004m.w4[is.na(sc4$birthM)]
sc4$birthM[is.na(sc4$birthM)] <- sc4$t70004m.w5[is.na(sc4$birthM)]
sc4$birthM[is.na(sc4$birthM)] <- sc4$t70004m.w6[is.na(sc4$birthM)]
starttime <- 12*(2010-1900)+9
cmcSt <- 12*(sc4$birthY-1900)+sc4$birthM
sc4$age <- (starttime-cmcSt)/12

# generation status (0 = no migration, 1 = migration, migration yes if gen.status < 3)
sc4$migration                       <- sc4$t400500_g1.w1
sc4$migration[is.na(sc4$migration)] <- sc4$t400500_g1.w2[is.na(sc4$migration)]
sc4$migration[is.na(sc4$migration)] <- sc4$t400500_g1.w3[is.na(sc4$migration)]
sc4$migration[is.na(sc4$migration)] <- sc4$t400500_g1.w4[is.na(sc4$migration)]
sc4$migration[is.na(sc4$migration)] <- sc4$t400500_g1.w5[is.na(sc4$migration)]
sc4$migration[is.na(sc4$migration)] <- sc4$t400500_g1.w6[is.na(sc4$migration)]
sc4$migration <- ifelse(is.na(sc4$migration), NA, ifelse(sc4$migration %in% 1:6, 1, 0))

# competencies
sc4$math1 <- sc4$mag9_sc1
sc4$math2 <- sc4$mag12_sc1u

# DGCF
sc4$dgcf1 <- recodeVar(sc4$dgci2101_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf2 <- recodeVar(sc4$dgci2102_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf3 <- recodeVar(sc4$dgci2103_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf4 <- recodeVar(sc4$dgci2104_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf5 <- recodeVar(sc4$dgci2201_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf6 <- recodeVar(sc4$dgci2202_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf7 <- recodeVar(sc4$dgci2203_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf8 <- recodeVar(sc4$dgci2204_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf9 <- recodeVar(sc4$dgci2301_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf10 <- recodeVar(sc4$dgci2302_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf11 <- recodeVar(sc4$dgci2303_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf12 <- recodeVar(sc4$dgci2304_sc4g9_c, c(0, 1), c(0, 1), default = NA)
sc4$dgcf <- rowSums(sc4[, paste0("dgcf", 1:12)], na.rm = TRUE)
table(sc4$dgcf, sc4$dgg9_sc3b, useNA = "always")
sc4$dgcf[is.na(sc4$dgg9_sc3b)] <- NA

# reliability: DGCF
ci.reliability(sc4[, paste0("dgcf", 1:12)], type = "categorical",
               interval.type = "none")

# mathematical self-concept
sc4$sc1 <- recodeVar(sc4$t66001a, 1:4, 1:4, default = NA) # I get good grades in math                       
sc4$sc2 <- recodeVar(sc4$t66001b, 1:4, 1:4, default = NA) # Math is one of my best subjects.
sc4$sc3 <- recodeVar(sc4$t66001c, 1:4, 1:4, default = NA) # I have always been good in Math
sc4$sc <- rowMeans(sc4[, paste0("sc", 1:3)], na.rm = TRUE)
sc4$sc[is.nan(sc4$sc)] <- NA

# reliability: sc
ci.reliability(sc4[, paste0("sc", 1:3)], interval.type = "none")

# school type (gen)
sc4$schooltype_w1 <- recodeVar(as.numeric(sc4$t723080_g1.w1), c(3:9, 11:17),
  c("HS", "RS", rep("REST",4), "GY", rep("REST",7)), default = NA)
sc4$schooltype_exp <- recodeVar(as.numeric(sc4$stratum_exp), c(1, 2, 3, 4:6), 
  c("GY", "RS", "HS", rep("REST",3)))   
sc4$schooltype_w1[is.na(sc4$schooltype_w1)] <- sc4$schooltype_exp[is.na(sc4$schooltype_w1)]  

# -------------------------------------------------------------------------
# ADD INFORMATION ON ASSESSMENT MODE (provided by NEPS method group)
# -------------------------------------------------------------------------

sc4 <- merge(sc4, dat6, by="ID_t", all.x=TRUE)
sc4$amode <- NA
sc4$amode[sc4$tstud_ak_7 %in% "A50"] <- 0 # "school context"
sc4$amode[!(sc4$tstud_ak_7 %in% "A50")] <- 1 # "individual field"

# --------------------------------------------------------------------
# SELECT SUBSAMPLE
# --------------------------------------------------------------------

# participants at wave 1
items <- c("ID_t", "ID_i",
           "math1", "math2",
           "dgcf", "sc", 
           "sex", "age", "migration",
           "schooltype_w1", "amode")
sc4 <- sc4[as.numeric(sc4$wave_w1) == 1, items]
rm(items)

# valid mathematical score at wave 1
sc4 <- subset(sc4, !is.na(math1))

# 1106 missings in dgcf variable: 7.6% of cases -> take them out
sc4 <- subset(sc4, !is.na(dgcf))   

# percentage of missing values on covariates (TRUE = missing) -> negligible
prop.table(table(is.na(sc4$dgcf)))
prop.table(table(is.na(sc4$sc)))
prop.table(table(is.na(sc4$sex)))
prop.table(table(is.na(sc4$age)))
prop.table(table(is.na(sc4$migration)))
prop.table(table(is.na(sc4$schooltype_w1)))
prop.table(table(is.na(sc4$amode)))

# --------------------------------------------------------------------
# HANDLING MISSING VALUES IN COVARIATES BY SIMPLE IMPUTATION
# --------------------------------------------------------------------
datSel <- sc4[,c("math1","dgcf", "sc", "sex", "age", "migration", "schooltype_w1", "amode")]
datSel$schooltype_w1 <- as.factor(datSel$schooltype_w1)
datSel$schooltype_w1  <- relevel(datSel$schooltype_w1, ref="GY")
md.pattern(datSel) # missing cases for migration & sc less than 1 percent: single imputation with MICE
datI <- mice(datSel, m = 1, seed = 235)
datIm <- complete(datI, action = 1)
sc4 <- cbind(sc4[, c("ID_t", "ID_i", "math2")],datIm)
md.pattern(sc4)

# --------------------------------------------------------------------
# DESCRIPTIVES
# --------------------------------------------------------------------

sc4$HS <- ifelse(sc4$schooltype_w1 %in% "HS", 1, 0)
sc4$RS <- ifelse(sc4$schooltype_w1 %in% "RS", 1, 0)
sc4$REST <-  ifelse(sc4$schooltype_w1 %in% "REST", 1, 0)

# descriptives
apply(sc4[, c("math1", "math2", "sc", "dgcf", "sex", "migration", "age", "amode","HS","RS", "REST")], 2,
      function(x) { round(c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)), 3) })

# correlations
round(cor(sc4[, c("math1", "math2", "sc", "dgcf",
                    "sex", "migration", "age", "amode", "HS", "RS", "REST")], use = "pairwise"), 3)

# --------------------------------------------------------------------
# ASSESS SELECTION BIAS
# --------------------------------------------------------------------
compNiv1 <- quantile(sc4$math1)
Smiss <- sc4[is.na(sc4$math2),]
table(cut(Smiss$math1,compNiv1))
# 62% lower than median compet. Wave 1

# --------------------------------------------------------------------
# LAST DATA MANIPULATIONS  
# --------------------------------------------------------------------

sc4$sc <- sc4$sc-mean(sc4$sc)             # z-standardizing
sc4$sc <- sc4$sc/sqrt(var(sc4$sc))
sc4$dgcf <- sc4$dgcf-mean(sc4$dgcf)
sc4$dgcf <- sc4$dgcf/sqrt(var(sc4$dgcf))
sc4$DROP <- ifelse(is.na(sc4$math2), 1,0)   # indicator for dropping out

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}       # make everything a numerical value
as.numeric.matrix <- function(mat){
 for(i in 1:ncol(mat)){
    if(is.factor(mat[,i]))
      mat[,i] <- as.numeric.factor(mat[,i])
 }
 return(as.matrix(mat))
}

sdat <- as.numeric.matrix(sc4[, c("ID_t","ID_i","dgcf","sc","sex","HS","RS","REST","migration","amode","age","math1","math2","DROP")])

# --------------------------------------------------------------------
# SAVE DATA
# --------------------------------------------------------------------

write.table(x=sdat, sep=",", col.names=FALSE, row.names=FALSE, na=".",
 file="Z:\\Projektgruppen_(08)\\StatMethoden_(p000033)\\SelectivityLinPanelModel\\ModelCodeVersion1\\syntax\\SupplementalMaterialPaper\\dataSC4MatheWIDE.csv")



