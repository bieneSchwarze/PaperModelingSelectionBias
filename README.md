# PaperModelingSelectionBias

This material belongs to the paper draft "Modeling Competence Development in the Presence of Selection Bias" by S. Zinn and T. Gnambs.

It contains the source code of statistical software used to conduct sensitivity analyses when dealing with missing data in longitudinal and multilevel modeling. Concretely, the competence development of German adolescents is modelled. 

The following models have been estimated: listwise deletion LWD, full information maximum likelihood FIML, weighted regression with inverse probabiltiy weights WE, mltivariate imputation via chained equations MI, Diggle-Kenward model DK, Wu-Carrol model WC, and Roy pattern mixture model PM.

The programs in this repository are named as follows:
Data preparation and description:	Prepare+describeData.r	(R file),
Listwise deletion model: LWD.r	(R file),
Full information maximum likelihood: 	FIML.inp	(Mplus file),
Nonresponse analysis, computation of intra class correlation, and inverse probability weighting: NonresponseAnalysis+ICC+WE.do	(Stata file),
Little test and multivariate imputation via chained equations: Little+MI.r	(R file),
Diggle-Kenward selection model: DK.inp	(Mplus file) ,
Wu-Carroll  selection model:  WC.inp	(Mplus),
Roy pattern mixture model: 	PM.inp	(Mplus).

Beware that using the data from the National Educational Panel Study (NEPS): Starting Cohort Grade 9, doi:10.5157/NEPS:SC4:9.0.0 as done in this study requires a contract of use with the Leibniz Institute for Educational Trajectories in Bamberg/Germany. 
