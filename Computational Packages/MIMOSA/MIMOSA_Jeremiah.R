library(MIMOSA)
library(dplyr)
library(data.table)

###Here is the updated script
setwd("/Users/tarynam/Desktop/Jeremiah/")
#I manually edited these column names in the csv file because they were horrible to format
DF<-read.csv("JKK Flow Data Absolute Counts 2019_AUG_19.csv")
names(DF)<-gsub("CD4_CD4.", "CD4_", names(DF))
names(DF)<-gsub("CD8_CD8.", "CD8_", names(DF))
a<-tstrsplit(DF$Sample, " ")[[5]]
b<-tstrsplit(a, "_")[[1]]
b<-gsub(",2f,", ".", b)
DF$condition<-b


#CD8 IFNg
bubblegum<-dplyr::select(DF, PATIENT.ID, condition, CD8_IFNgP, CD8_Total) %>%
    dplyr::mutate(Neg = (CD8_Total - CD8_IFNgP)) %>%
    dplyr::select(-CD8_Total)

doozie<-ConstructMIMOSAExpressionSet(bubblegum,
                                     reference=condition %in% 'uns',
                                     measure.columns=c('Neg','CD8_IFNgP'),
                                     other.annotations=c('condition','PATIENT.ID'),
                                     default.cast.formula=component~PATIENT.ID+condition,
                                     .variables=.(PATIENT.ID))

banana <- MIMOSA(Neg+CD8_IFNgP~PATIENT.ID|condition,
                 data=doozie,
                 method="mcmc")

probabilities<-MIMOSA::getZ(banana)
counts<-MIMOSA::countsTable(banana)
fdr<-MIMOSA::fdr(banana)
MIMOSA.Results_ifng<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary_ifng<-MIMOSA::getW(banana)
MIMOSA.Results_ifng$measure<-"ifng"
MIMOSA.Results_ifng$Sample<-rownames(MIMOSA.Results_ifng)
names(MIMOSA.Results_ifng)<-gsub("CD8_IFNgP","Pos", names(MIMOSA.Results_ifng))
MIMOSA.Results_ifng$Sample<-gsub(".ESAT6", "", MIMOSA.Results_ifng$Sample)
#comapre to cheryl
cheryl<-read.csv("JKK_CD8_IFNg_MIMOSA_2020_MAR_26.csv")
cheryl$Sample<-gsub(".ESAT6", "", cheryl$Sample)
test<-merge(cheryl, MIMOSA.Results_ifng, by="Sample")
test$Neg.x-test$Neg.y
test$IFNg-test$Pos
test$IFNg_REF-test$Pos_REF
test$Neg_REF.x-test$Neg_REF.y
test$pr.nonresponse_dif<-abs(100*(test$Pr.Nonresponse.x-test$Pr.Nonresponse.y))
test$pr.response_dif <- abs(100*(test$Pr.response.x-test$Pr.response.y))
test$fdr_dif<-abs(100*(test$fdr.x-test$fdr.y))
test$cheryl<-0
test$cheryl[test$Pr.response.x>.9 & test$fdr.x<0.03]<-1
test$taryn<-0
test$taryn[test$Pr.response.y>.9 & test$fdr.y<0.03]<-1
test$agree<-"disagree"
test$agree[test$cheryl==1 & test$taryn==1]<-"agree"
test$agree[test$cheryl==0 & test$taryn==0]<-"agree"

a<-ggplot(test, aes(x=pr.response_dif, fill=agree))+
    geom_histogram()+facet_grid(~agree)+
    labs(title="Difference in Probability of Response (%)", x="")
b<-ggplot(test, aes(x=fdr_dif, fill=agree))+
    geom_histogram()+facet_grid(~agree)+
    labs(title="Difference in FDR (%)", x="")

c<-ggplot(test, aes(x=Pr.response.x, y=Pr.response.y, color=agree))+
    geom_point()+geom_smooth()+facet_grid(~agree)+
    labs(x="Cheryl's Values", y="Taryn's values",
         title="Probability of Response")
d<-ggplot(test, aes(x=fdr.x, y=fdr.y, color=agree))+
    geom_point()+geom_smooth()+facet_grid(~agree)+
    labs(x="Cheryl's Values", y="Taryn's values",
         title="FDR")

cowplot::plot_grid(a, b, c, d, nrow=2)



#CD8 CD107
bubblegum<-dplyr::select(DF, PATIENT.ID, condition, CD8_CD107aP, CD8_Total) %>%
    dplyr::mutate(Neg = (CD8_Total - CD8_CD107aP)) %>%
    dplyr::select(-CD8_Total)

doozie<-ConstructMIMOSAExpressionSet(bubblegum,
                                     reference=condition %in% 'uns',
                                     measure.columns=c('Neg','CD8_CD107aP'),
                                     other.annotations=c('condition','PATIENT.ID'),
                                     default.cast.formula=component~PATIENT.ID+condition,
                                     .variables=.(PATIENT.ID))

banana <- MIMOSA(Neg+CD8_CD107aP~PATIENT.ID|condition,
                 data=doozie,
                 method="mcmc")

probabilities<-MIMOSA::getZ(banana)
counts<-MIMOSA::countsTable(banana)
fdr<-MIMOSA::fdr(banana)
MIMOSA.Results_107a<-data.frame(cbind(counts,probabilities,fdr))
MIMOSA.summary_107a<-MIMOSA::getW(banana)
MIMOSA.Results_107a$measure<-"107a"
MIMOSA.Results_107a$Sample<-rownames(MIMOSA.Results_107a)
names(MIMOSA.Results_107a)<-gsub("CD8_CD107aP","Pos", names(MIMOSA.Results_107a))

test<-rbindlist(list(MIMOSA.Results_ifng, MIMOSA.Results_107a))
write.csv(test, "CD8results.csv")
