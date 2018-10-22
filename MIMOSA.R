biocLite("devtools")
library(devtools) #load it
install_github("MIMOSA","RGLab",branch="master") #Need R 3.0.0, and a bunch of dependencies. The install will fail with various error messages until you install those dependencies.
biocLite("MIMOSA")
library(MIMOSA)

DF<-read.csv("/Applications/Old Computer/Day Lab/Data_Flow/Clean/All Counts.csv")
DF<-as.data.table(DF) 
DF$Neg<-DF[,21]
DF<-DF[DF$Stim!="PMA"]

#All cytokine
data<-DF[,c(2,5,26,22)]
E<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim%in%'UN',
                                measure.columns=c('Neg','Cyt.Count'),
                                other.annotations=c('Stim','ID'),
                                default.cast.formula=component~ID+Stim,
                                .variables=.(ID))
test <- MIMOSA(Neg+Cyt.Count~ID|Stim,
               data=E)
probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results<-cbind(counts,probabilities,fdr)
MIMOSA.summary<-MIMOSA::getW(test)
write.csv(MIMOSA.Results,"/Applications/Old Computer/Day Lab/Data_Flow/Clean/MIMOSA-Results.csv")

#TH1 
data1<-DF[,c(2,5,26,23)]
E1<-ConstructMIMOSAExpressionSet(data1,
                                reference=Stim%in%'UN',
                                measure.columns=c('Neg','TH1.Count'),
                                other.annotations=c('Stim','ID'),
                                default.cast.formula=component~ID+Stim,
                                .variables=.(ID))

test1 <- MIMOSA(Neg+TH1.Count~ID|Stim,
                data=E1)
probabilities1<-MIMOSA::getZ(test1)
counts1<-MIMOSA::countsTable(test1)
fdr1<-MIMOSA::fdr(test1)
MIMOSA.Results1<-cbind(counts1,probabilities1,fdr1)
MIMOSA.summary1<-MIMOSA::getW(test1)
write.csv(MIMOSA.Results1,"/Applications/Old Computer/Day Lab/Data_Flow/Clean/MIMOSA-Results1.csv")


#TH2 
data2<-DF[,c(2,5,26,25)]
E2<-ConstructMIMOSAExpressionSet(data2,
                                 reference=Stim%in%'UN',
                                 measure.columns=c('Neg','TH2.Count'),
                                 other.annotations=c('Stim','ID'),
                                 default.cast.formula=component~ID+Stim,
                                 .variables=.(ID))

test2 <- MIMOSA(Neg+TH2.Count~ID|Stim,
                data=E2)
probabilities2<-MIMOSA::getZ(test2)
counts2<-MIMOSA::countsTable(test2)
fdr2<-MIMOSA::fdr(test2)
MIMOSA.Results2<-cbind(counts2,probabilities2,fdr2)
MIMOSA.summary2<-MIMOSA::getW(test)
write.csv(MIMOSA.Results2,"/Applications/Old Computer/Day Lab/Data_Flow/Clean/MIMOSA-Results2.csv")


#TH1/2 
data1.2<-DF[,c(2,5,26,24)]
E1.2<-ConstructMIMOSAExpressionSet(data1.2,
                                 reference=Stim%in%'UN',
                                 measure.columns=c('Neg','TH1.2.Count'),
                                 other.annotations=c('Stim','ID'),
                                 default.cast.formula=component~ID+Stim,
                                 .variables=.(ID))

test1.2 <- MIMOSA(Neg+TH1.2.Count~ID|Stim,
                data=E1.2)
probabilities1.2<-MIMOSA::getZ(test1.2)
counts1.2<-MIMOSA::countsTable(test1.2)
fdr1.2<-MIMOSA::fdr(test1.2)
MIMOSA.Results1.2<-cbind(counts1.2,probabilities1.2,fdr1.2)
MIMOSA.summary1.2<-MIMOSA::getW(test1.2)
write.csv(MIMOSA.Results1.2,"/Applications/Old Computer/Day Lab/Data_Flow/Clean/MIMOSA-Results12.csv")


#real examples
E2<-ConstructMIMOSAExpressionSet(ICS,
                                 reference=ANTIGEN%in%'negctrl',
                                 measure.columns=c('CYTNUM','NSUB'),
                                 other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
                                 default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
                                 .variables=.(TCELLSUBSET,CYTOKINE,UID),
                                 featureCols=1,ref.append.replace='_REF')


result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
               data=E2, 
               method='EM',
               subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
               ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')


# add this in
test.TH2$Positive<-NA
for (i in 1:length(test.TH2$Positive)){
  if(test.TH2$fdr[i]<0.05){
    test.TH2$Positive[i] <- "Pos"}

  else{
    test.TH2$Positive[i] <- "Neg"} 
}

###(.fitMCMC ...EXPRATE ... currently set to 10E-4 ... set it yourself to 1-10% ... .01 to .05)

##set method to MCMC
