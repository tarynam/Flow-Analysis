install.packages("remotes") #install devtools package
library(remotes) #load it
install_github("RGLab/MIMOSA",ref="trunk") 
#Need R 3.0.0, and a bunch of dependencies.
#The install will fail with various error messages until you install those dependencies.
library(MIMOSA)
library(dplyr)
library(data.table)

###Here is the updated script
setwd("/Users/tarynam/Desktop/")
DF<-read.csv("JKK Flow Data Absolute Counts 2019_AUG_19.csv")
#show data table --> column names
DF$cd40 <- DF$CD4.40.p_107.p_G.p + DF$CD4.40.p_107.p_G.n + 
    DF$CD4.40.p_107.n_G.p + DF$CD4.40.p_107.n_G.n
DF$total<- DF$CD4.40.p_107.p_G.p + DF$CD4.40.p_107.p_G.n + 
    DF$CD4.40.p_107.n_G.p + DF$CD4.40.p_107.n_G.n + 
    DF$CD4.40.n_107.p_G.p + DF$CD4.40.n_107.p_G.n + 
    DF$CD4.40.n_107.n_G.p + DF$CD4.40.n_107.n_G.n
##easier for me this way
CD4<-mutate(DF,
            #These all follow the same format
            #IFNg is the new variable. it gets assign the row sum of all the columns that have a G in the name
            cd40=  rowSums(DF[,grep("CD4.40.p_107.._G..", names(DF))]),
            cd107=  rowSums(DF[,grep("CD4.40.._107.p_G..", names(DF))]),
            ifng=  rowSums(DF[,grep("CD4.40.._107.._G.p", names(DF))]),
            total=   rowSums(DF[,grep("CD4.40.._107.._G..", names(DF))])
)
CD4$condition <- tstrsplit(CD4$TUBE.NAME, " ")[[3]]
CD4$condition<-gsub("CFP10_ESAT6", "CFP10/ESAT6", CD4$condition)
names(CD4)

#IFNg
bubblegum<-dplyr::select(CD4, PATIENT.ID, condition, ifng, total) %>%
    dplyr::mutate(Neg = (total - ifng)) %>%
    dplyr::select(-total)

doozie<-ConstructMIMOSAExpressionSet(bubblegum,
                                reference=condition %in% 'uns',
                                measure.columns=c('Neg','ifng'),
                                other.annotations=c('condition','PATIENT.ID'),
                                default.cast.formula=component~PATIENT.ID+condition,
                                .variables=.(PATIENT.ID))

banana <- MIMOSA(Neg+ifng~PATIENT.ID|condition,
               data=doozie,
               method="mcmc")

probabilities<-MIMOSA::getZ(banana)
counts<-MIMOSA::countsTable(banana)
fdr<-MIMOSA::fdr(banana)
MIMOSA.Results_ifng<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary_ifng<-MIMOSA::getW(banana)


#cd40
bubblegum<-dplyr::select(CD4, PATIENT.ID, condition, cd40, total) %>%
    dplyr::mutate(Neg = (total - cd40)) %>%
    dplyr::select(-total)

doozie<-ConstructMIMOSAExpressionSet(bubblegum,
                                     reference=condition %in% 'uns',
                                     measure.columns=c('Neg','cd40'),
                                     other.annotations=c('condition','PATIENT.ID'),
                                     default.cast.formula=component~PATIENT.ID+condition,
                                     .variables=.(PATIENT.ID))

banana <- MIMOSA(Neg+cd40~PATIENT.ID|condition,
                 data=doozie,
                 method="mcmc")

probabilities<-MIMOSA::getZ(banana)
counts<-MIMOSA::countsTable(banana)
fdr<-MIMOSA::fdr(banana)
MIMOSA.Results_cd40<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary_cd40<-MIMOSA::getW(banana)

#cd107
bubblegum<-dplyr::select(CD4, PATIENT.ID, condition, cd107, total) %>%
    dplyr::mutate(Neg = (total - cd107)) %>%
    dplyr::select(-total)

doozie<-ConstructMIMOSAExpressionSet(bubblegum,
                                     reference=condition %in% 'uns',
                                     measure.columns=c('Neg','cd107'),
                                     other.annotations=c('condition','PATIENT.ID'),
                                     default.cast.formula=component~PATIENT.ID+condition,
                                     .variables=.(PATIENT.ID))

banana <- MIMOSA(Neg+cd107~PATIENT.ID|condition,
                 data=doozie,
                 method="mcmc")

probabilities<-MIMOSA::getZ(banana)
counts<-MIMOSA::countsTable(banana)
fdr<-MIMOSA::fdr(banana)
MIMOSA.Results_cd107<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary_cd107<-MIMOSA::getW(banana)


MIMOSA.Results_cd40$param<-"cd40"
MIMOSA.Results_cd40$sample<-rownames(MIMOSA.Results_cd40)

MIMOSA.Results_cd107$param<-"cd107"
MIMOSA.Results_cd107$sample<-rownames(MIMOSA.Results_cd107)

MIMOSA.Results_ifng$param<-"ifng"
MIMOSA.Results_ifng$sample<-rownames(MIMOSA.Results_ifng)

total<-rbindlist(list(MIMOSA.Results_ifng, MIMOSA.Results_cd40, MIMOSA.Results_cd107))
total$sample<-gsub("_Treatment", "", total$sample)
write.csv(total, "mimosa_forcheryl.csv")
