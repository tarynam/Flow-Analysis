install.packages("devtools") #install devtools package
library(devtools) #load it
install_github("MIMOSA","RGLab",branch="master") #Need R 3.0.0, and a bunch of dependencies. The install will fail with various error messages until you install those dependencies.
library(MIMOSA)

setwd("/Applications/Old Computer/Day Lab/Flow-Data/clean/")
DF<-read.csv("CD4_Counts_Clean.csv")

#All cytokine
data<-dplyr::select(DF, Sample, Donor, TB, SM, Stim, Cyt, Total) %>%
    dplyr::mutate(Neg = (Total - Cyt)) %>%
    dplyr::select(-Total)

E<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim %in% 'UN',
                                measure.columns=c('Neg','Cyt'),
                                other.annotations=c('Stim','Sample', 'Donor', 'TB', 'SM'),
                                default.cast.formula=component~Donor+Stim,
                                .variables=.(Donor))
test <- MIMOSA(Neg+Cyt~Donor|Stim,
               data=E,
               method="mcmc")

probabilities<-MIMOSA::getZ(test)
counts<-MIMOSA::countsTable(test)
fdr<-MIMOSA::fdr(test)
MIMOSA.Results<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary<-MIMOSA::getW(test)

#TH1 
data<-dplyr::select(DF, Sample, Donor, TB, SM, Stim, TH1, Total) %>%
    dplyr::mutate(Neg = (Total - TH1)) %>%
    dplyr::select(-Total)

E1<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim %in% 'UN',
                                measure.columns=c('Neg','TH1'),
                                other.annotations=c('Stim','Sample', 'Donor', 'TB', 'SM'),
                                default.cast.formula=component~Donor+Stim,
                                .variables=.(Donor))

test1 <- MIMOSA(Neg+TH1~Donor|Stim,
                data=E1,
                method="mcmc")

probabilities1<-MIMOSA::getZ(test1)
counts1<-MIMOSA::countsTable(test1)
fdr1<-MIMOSA::fdr(test1)
MIMOSA.Results1<-data.frame(cbind(counts1,probabilities1,fdr1))
MIMOSA.summary1<-MIMOSA::getW(test1)


#TH2 
data<-dplyr::select(DF, Sample, Donor, TB, SM, Stim, TH2, Total) %>%
    dplyr::mutate(Neg = (Total - TH2)) %>%
    dplyr::select(-Total)

E2<-ConstructMIMOSAExpressionSet(data,
                                 reference=Stim %in% 'UN',
                                 measure.columns=c('Neg','TH2'),
                                 other.annotations=c('Stim','Sample', 'Donor', 'TB', 'SM'),
                                 default.cast.formula=component~Donor+Stim,
                                 .variables=.(Donor))

test2 <- MIMOSA(Neg+TH2~Donor|Stim,
                data=E2,
                method="mcmc")

probabilities2<-MIMOSA::getZ(test2)
counts2<-MIMOSA::countsTable(test2)
fdr2<-MIMOSA::fdr(test2)
MIMOSA.Results2<-data.frame(cbind(counts2,probabilities2,fdr2))
MIMOSA.summary2<-MIMOSA::getW(test)

#TH1/2 
data<-dplyr::select(DF, Sample, Donor, TB, SM, Stim, TH1.2, Total) %>%
    dplyr::mutate(Neg = (Total - TH1.2)) %>%
    dplyr::select(-Total)

E1.2<-ConstructMIMOSAExpressionSet(data,
                                 reference=Stim %in% 'UN',
                                 measure.columns=c('Neg','TH1.2'),
                                 other.annotations=c('Stim','Sample', 'Donor', 'TB', 'SM'),
                                 default.cast.formula=component~Donor+Stim,
                                 .variables=.(Donor))

test1.2 <- MIMOSA(Neg+TH1.2~Donor|Stim,
                data=E1.2,
                method="mcmc")
probabilities1.2<-MIMOSA::getZ(test1.2)
counts1.2<-MIMOSA::countsTable(test1.2)
fdr1.2<-MIMOSA::fdr(test1.2)
MIMOSA.Results1.2<-data.frame(cbind(counts1.2,probabilities1.2,fdr1.2))
MIMOSA.summary1.2<-MIMOSA::getW(test1.2)
write.csv(MIMOSA.Results1.2,"/Applications/Old Computer/Day Lab/Data_Flow/Clean/MIMOSA-Results12.csv")

names(MIMOSA.Results)<-paste("Cyt", names(MIMOSA.Results), sep="_")
names(MIMOSA.Results1)<-paste("TH1", names(MIMOSA.Results1), sep="_")
names(MIMOSA.Results2)<-paste("TH2", names(MIMOSA.Results2), sep="_")
names(MIMOSA.Results1.2)<-paste("TH1.2", names(MIMOSA.Results1.2), sep="_")

MIMOSA.Results$Sample <- row.names(MIMOSA.Results)
MIMOSA.Results1$Sample <- row.names(MIMOSA.Results1)
MIMOSA.Results2$Sample <- row.names(MIMOSA.Results2)
MIMOSA.Results1.2$Sample <- row.names(MIMOSA.Results1.2)

dfs<-list(MIMOSA.Results, MIMOSA.Results1, MIMOSA.Results2, MIMOSA.Results1.2)
mimosa<-join_all(dfs, "Sample")

mimosa$Cyt_include<-0
for(i in 1:nrow(mimosa)){
    if(mimosa$Cyt_Pr.response[i]>0.7 & mimosa$Cyt_fdr[i]<0.03){
        mimosa$Cyt_include[i]<-1
        print(mimosa$Cyt_include[i])}
}

mimosa$TH1_include<-0
for(i in 1:nrow(mimosa)){
    if(mimosa$TH1_Pr.response[i]>0.7 & mimosa$TH1_fdr[i]<0.03){
        mimosa$TH1_include[i]<-1
        print(mimosa$TH1_include[i])}
}

mimosa$TH2_include<-0
for(i in 1:nrow(mimosa)){
    if(mimosa$TH2_Pr.response[i]>0.7 & mimosa$TH2_fdr[i]<0.03){
        mimosa$TH2_include[i]<-1
        print(mimosa$TH2_include[i])}
}

mimosa$TH1.2_include<-0
for(i in 1:nrow(mimosa)){
    if(mimosa$TH1.2_Pr.response[i]>0.7 & mimosa$TH1.2_fdr[i]<0.03){
        mimosa$TH1.2_include[i]<-1
        print(mimosa$TH1.2_include[i])}
}

mimosa$score<-(mimosa$Cyt_include + mimosa$TH1_include + mimosa$TH2_include + mimosa$TH1.2_include)


mimosa$Sample<-gsub("_Treatment", "", mimosa$Sample)
split<-tstrsplit(mimosa$Sample, "_")
mimosa$Donor<-split[[1]]
mimosa$Stim<-split[[2]]
mimosa<-dplyr::select(mimosa, Sample, Donor, Stim, everything())
write.csv(mimosa, "CD4_MIMOSA_Complete.csv")



