setwd("/Users/oldmac/Desktop/Day Lab/Data")
NK_ICS<-read.csv("NK_ICS.csv")
OG<-read.csv("EVC_OG.csv")
IFN<-read.csv("ELISA_IFN.csv")
IL13<-read.csv("ELISA_IL13.csv")

##ELISA###
IFN$Cytokine<-as.character("IFNg")
IL13$Cytokine<-as.character("13")
setnames(IFN,old="Dilution.Factor",new="Dilution")
ELISA<-rbind(IFN,IL13)
splStudyID = strsplit(as.character(ELISA$Sample)," ");
ELISA$Donor = as.factor(sapply(splStudyID,"[",1));
ELISA$Antigen <- as.factor(sapply(splStudyID, "[", 2));
ELISA$Result[ELISA$Result>500]<-500
ELISA$Adjusted<-(ELISA$Result*ELISA$Dilution)

## combine and melt
melted<-melt(ELISA, id.vars=c("Donor", "Antigen", "Cytokine"))
#for some reason it stores all the values as character so coerce them to numeric to run stats
melted$value<-as.numeric(as.character(melted$value))
stats<-as.data.frame(
  ddply(melted, c("Antigen", "Donor", "Cytokine", "variable"), summarise,
        mean = mean(value), sd = sd(value),
        sem = sd(value)/sqrt(length(value))))

# re-separate data sets
ELISA.Mean<-subset(stats,stats$variable=="Adjusted")
ELISA.Mean<-ELISA.Mean[!("Antigen"=="Pep")]


ELISA.Mean[,"Group"]<-NA
ELISA.Mean$Group[ELISA.Mean$Donor=="HD235"]<-as.character("TH2")
ELISA.Mean$Group[ELISA.Mean$Donor=="HD355"]<-as.character("HD")
ELISA.Mean$Group[ELISA.Mean$Donor=="HD203"]<-as.character("HD")
ELISA.Mean$Group[ELISA.Mean$Donor=="PS16-1001"]<-as.character("SA")
ELISA.Mean$Group[ELISA.Mean$Donor=="PS16-1002"]<-as.character("SA")
ELISA.Mean$Group[is.na(ELISA.Mean$Group)]<-as.character("SM")

ELISA.Mean<-ELISA.Mean[!(Donor=="NK2020")]

IFNg<-subset(ELISA.Mean,ELISA.Mean$Cytokine=="IFNg")
IFNg<-subset(ELISA.Mean,ELISA.Mean$Group=="SM")
ggplot(IFNg, aes(x=Antigen, y=mean, color=Donor)) +
  geom_point(shape=16, size=5)+
  xlab("Antigen") +
  ylab("IFNg pg/mL") +
  ggtitle("ELISA") +
  theme_bw()+
  theme(text = element_text(size = 20))

IL13<-subset(ELISA.Mean,ELISA.Mean$Cytokine=="13")
IL13<-subset(ELISA.Mean, ELISA.Mean$Group=="SM")
IL13<-IL13[!(Antigen=="Pep")]

ggplot(IL13, aes(x=Antigen, y=mean, color=Donor)) +
  geom_point(shape=16, size=5)+
  xlab("Antigen") +
  ylab("IL13 pg/mL") +
  ggtitle("ELISA") +
  theme_bw()+
  theme(text = element_text(size = 20))

##OG##
splStudyID = strsplit(as.character(OG$Sample),":");
OG$Donor = as.factor(sapply(splStudyID,"[",2));

splStudyID = strsplit(as.character(OG$Donor),"_");
OG$Donor = as.factor(sapply(splStudyID,"[",1));
OG$Antigen = as.factor(sapply(splStudyID, "[",2))
OG$Stim = factor(OG$Antigen,levels=c("UN","PPD","WCL","PEP","SEB"))

OG[,"GATA3"]<-NA
for (i in 1:length(OG$GATA3)){
  OG$GATA3[i] <- ifelse((OG$Stim[i]=="UN"), as.numeric(OG$CD4.GATA3[i]),
                                as.numeric(OG$OG.GATA3[i]))
}
OG$GATA3<-as.numeric(OG$GATA3)

OG[,"TH1.2_TF"]<-NA
for (i in 1:length(OG$TH1.2_TF)){
  OG$TH1.2_TF[i] <- ifelse((OG$Stim[i]=="UN"), as.numeric(OG$CD4.TH1.TH2_TF[i]),
                        as.numeric(OG$OG.TH1.TH2_TF[i]))
}
OG$TH1.2_TF<-as.numeric(OG$TH1.2_TF)

OG[,"TH1.2_SM"]<-NA
for (i in 1:length(OG$TH1.2_SM)){
  OG$TH1.2_SM[i] <- ifelse((OG$Stim[i]=="UN"), as.numeric(OG$CD4.TH1.TH2_SM[i]),
                           as.numeric(OG$OG.TH1.TH2_SM[i]))
}
OG$TH1.2_SM<-as.numeric(OG$TH1.2_SM)

OG[,"Tbet"]<-NA
for (i in 1:length(OG$Tbet)){
  OG$Tbet[i] <- ifelse((OG$Stim[i]=="UN"), as.numeric(OG$CD4.Tbet[i]),
                           as.numeric(OG$OG.Tbet[i]))
}
OG$Tbet<-as.numeric(OG$Tbet)

####PLOTS
dodge <- position_dodge(width=0.5)  

#Proliferation
ggplot(OG, aes(x=Stim, y=CD4.OG, color=Donor)) +
  geom_point(shape=16, size=10, position=dodge) + 
  xlab("Stimulus") +
  ylab("Frequency of OG Low Cells") +
  ggtitle("Oregon Green Proliferation") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

#GATA3
ggplot(OG, aes(x=Stim, y=GATA3, color=Donor)) +
  geom_point(shape=16, size=10, position=dodge)+ 
  xlab("Stimulus") +
  ylab("Frequency of GATA3+ Cells") +
  ggtitle("GATA3 Expression in CD4+ T cells") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

#TF
ggplot(OG, aes(x=Stim, y=TH1.2_TF, color=Donor)) +
  geom_point(shape=16, size=10, position=dodge)+ 
  xlab("Stimulus") +
  ylab("Frequency of T-bet+GATA3+ Cells") +
  ggtitle("Transcription Factor Co-expression") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

#SM
ggplot(OG, aes(x=Stim, y=TH1.2_SM, color=Donor)) +
  geom_point(shape=16, size=10, position=dodge)+ 
  xlab("Stimulus") +
  ylab("Frequency of CCR4+CXCR3+ Cells") +
  ggtitle("Surface Marker Co-expression") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

#Tbet
ggplot(OG, aes(x=Stim, y=Tbet, color=Donor)) +
  geom_point(shape=16, size=10, position=dodge)+ 
  xlab("Stimulus") +
  ylab("Frequency of T-bet+ Cells") +
  ggtitle("T-bet Expression in CD4+ T Cells") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")


