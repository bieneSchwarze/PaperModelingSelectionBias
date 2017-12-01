********************************************************************************
********************************************************************************
*** SC4 - math competence
***
*** Analyses : Linear Model with Correction for Selection Bias
***
*** Method:  Weighting 
***
*** project Modeling longitudinal dropout
*** @author Sabine Zinn 
***
*** 16.08.2017 
***
********************************************************************************
********************************************************************************

import delimited "Z:\Projektgruppen_(08)\StatMethoden_(p000033)\SelectivityLinPanelModel\ModelCodeVersion1\syntax\SupplementalMaterialPaper\dataSC4MatheWIDE.csv", clear

rename v1 ID_t
rename v2 ID_i
rename v3 dgcf
rename v4 sc
rename v5 sex
rename v6 HS
rename v7 RS
rename v8 REST
rename v9 mig
rename v10 amode
rename v11 age
rename v12 math1
rename v13 math2
rename v14 dr

recode dr (0 = 1) (1 = 0), gen(tn)

// nonresponse analysis
gen sqmath = math1^2
xtmelogit dr math1 sqmath amode dgcf sc sex mig age HS RS REST || ID_i:

// create weights: only fixed effects
xtmelogit tn math1 sqmath amode dgcf sc sex mig age HS RS REST, noconstant || ID_i:

predict ree*, reffects // obtain the random effects
// only for 85% of the schools

predict xb, xb

gen muE = 1 /(1+exp(-1*(xb)))

gen wE = 1/muE

graph box wE tn

hist wE

// listwise deletion
//mixed math i.time c.sc c.dgcf c.sc#i.time c.dgcf#i.time || ID_i: || ID_t:, reml

// weighting
keep  if tn==1

hist wE

// long format
gen math2F = real(math2)
drop math2
rename math2F math2

reshape long math, i(ID_t)

recode _j (1 = 0) (2 = 1), gen(time)

// ICC
mixed math || ID_i:, || ID_t:, 

estat icc

mixed math i.time c.sc c.dgcf c.sc#i.time c.dgcf#i.time [pweight=wEl] || ID_i:, || ID_t:, 

estat ic
