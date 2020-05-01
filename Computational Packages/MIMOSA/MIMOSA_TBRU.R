install.packages("remotes") #install devtools package
library(remotes) #load it
install_github("RGLab/MIMOSA",ref="trunk") 
#Need R 3.0.0, and a bunch of dependencies.
#The install will fail with various error messages until you install those dependencies.
library(MIMOSA)
library(dplyr)
library(data.table)

###Here is the updated script
setwd("/Users/tarynam/Desktop/MIMOSA_docs_for_troubleshooting")
DF<-read.csv("TBRU Project 1_Counts_FOR COMPASS (1).csv")
#CD4.p is the total CD4 count
#CD4.p_IFNg.p is the IFNg postive count
#Sample_ID is the ID v

#IFNg
bubblegum<-dplyr::select(DF, Sample_ID, TimePoint, Stim, CD4.p, CD4.p_IFNg.p) %>%
    dplyr::mutate(Neg = (CD4.p - CD4.p_IFNg.p)) %>%
    dplyr::select(-CD4.p)

E<-ConstructMIMOSAExpressionSet(bubblegum,
                                reference=Stim %in% 'US',
                                measure.columns=c('Neg','CD4.p_IFNg.p'),
                                other.annotations=c('Stim','Sample_ID', 'TimePoint'),
                                default.cast.formula=component~Sample_ID+Stim+TimePoint,
                                .variables=.(Sample_ID, TimePoint))

results <- MIMOSA(Neg+CD4.p_IFNg.p~Sample_ID+TimePoint|Stim,
               data=E,
               method="mcmc",
               ref= RefTreat%in%'Reference')

probabilities<-MIMOSA::getZ(results)
counts<-MIMOSA::countsTable(results)
fdr<-MIMOSA::fdr(results)
MIMOSA.Results_ifng<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary_ifng<-MIMOSA::getW(results)


