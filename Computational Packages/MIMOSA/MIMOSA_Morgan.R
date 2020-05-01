biocLite("devtools")
library(devtools) #load it
install_github("MIMOSA","RGLab",branch="master") #Need R 3.0.0, and a bunch of dependencies. The install will fail with various error messages until you install those dependencies.
biocLite("MIMOSA")
library(MIMOSA)

DF<-read.csv("/Users/tarynam/Desktop/QFT-_HIV+_donor_file/CD8_TB_HIVexport_180824.csv")
DF<-as.data.table(DF) 


#GT
DF$Neg<-DF$COUNT-DF$GT_Count
DF<-DF[DF$CONDITION!="SEB"]
data<-dplyr::select(DF, DONOR, CONDITION, DISEASE, Neg, GT_Count)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=CONDITION%in%'unstim',
                                measure.columns=c('Neg','GT_Count'),
                                other.annotations=c('CONDITION','DONOR'),
                                default.cast.formula=component~DONOR+CONDITION,
                                .variables=.(DONOR))
test <- MIMOSA(Neg+GT_Count~DONOR|CONDITION,
               data=E, method="mcmc")
probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results1<-cbind(counts,probabilities,fdr)
MIMOSA.summary<-MIMOSA::getW(test)

#TNF-any
DF$Neg<-DF$COUNT-DF$TNFa_Count
DF<-DF[DF$CONDITION!="SEB"]
data<-dplyr::select(DF, DONOR, CONDITION, DISEASE, Neg, TNFa_Count)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=CONDITION%in%'unstim',
                                measure.columns=c('Neg','TNFa_Count'),
                                other.annotations=c('CONDITION','DONOR'),
                                default.cast.formula=component~DONOR+CONDITION,
                                .variables=.(DONOR))
test <- MIMOSA(Neg+TNFa_Count~DONOR|CONDITION,
               data=E, method="mcmc")
probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results2<-cbind(counts,probabilities,fdr)
MIMOSA.summary<-MIMOSA::getW(test)

#IFN- any
DF$Neg<-DF$COUNT-DF$IFNg_Count
DF<-DF[DF$CONDITION!="SEB"]
data<-dplyr::select(DF, DONOR, CONDITION, DISEASE, Neg, IFNg_Count)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=CONDITION%in%'unstim',
                                measure.columns=c('Neg','IFNg_Count'),
                                other.annotations=c('CONDITION','DONOR'),
                                default.cast.formula=component~DONOR+CONDITION,
                                .variables=.(DONOR))
test <- MIMOSA(Neg+IFNg_Count~DONOR|CONDITION,
               data=E, method="mcmc")
probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results3<-cbind(counts,probabilities,fdr)
MIMOSA.summary<-MIMOSA::getW(test)

#TNF-single
DF$TNF<-DF$TNFa_Count-DF$GT_Count
DF$Neg<-DF$COUNT-DF$TNF
DF<-DF[DF$CONDITION!="SEB"]
data<-dplyr::select(DF, DONOR, CONDITION, DISEASE, Neg, TNF)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=CONDITION%in%'unstim',
                                measure.columns=c('Neg','TNF'),
                                other.annotations=c('CONDITION','DONOR'),
                                default.cast.formula=component~DONOR+CONDITION,
                                .variables=.(DONOR))
test <- MIMOSA(Neg+TNF~DONOR|CONDITION,
               data=E, method="mcmc")
probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results4<-cbind(counts,probabilities,fdr)
MIMOSA.summary<-MIMOSA::getW(test)

#IFN- single
DF$IFN<-DF$IFNg_Count-DF$GT_Count
DF$Neg<-DF$COUNT-DF$IFN
DF<-DF[DF$CONDITION!="SEB"]
data<-dplyr::select(DF, DONOR, CONDITION, DISEASE, Neg, IFN)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=CONDITION%in%'unstim',
                                measure.columns=c('Neg','IFN'),
                                other.annotations=c('CONDITION','DONOR'),
                                default.cast.formula=component~DONOR+CONDITION,
                                .variables=.(DONOR))
test <- MIMOSA(Neg+IFN~DONOR|CONDITION,
               data=E, method="mcmc")
probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results5<-cbind(counts,probabilities,fdr)
MIMOSA.summary<-MIMOSA::getW(test)

