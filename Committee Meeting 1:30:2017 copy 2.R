setwd("/Users/oldmac/Desktop/Day Lab/Data")

######FIX PERM
FP<- read.csv("Fix Perm.csv")
splStudyID = strsplit(as.character(FP$Sample),":");
FP$Sample = as.factor(sapply(splStudyID,"[",2));
splStudyID = strsplit(as.character(FP$Sample),"_");
FP$Donor = as.factor(sapply(splStudyID,"[",1));
FP$Antigen = as.factor(sapply(splStudyID,"[",2));
FP$FP = as.factor(sapply(splStudyID,"[",3));
FP$Stim = factor(FP$Antigen,levels=c("UN", "WCL","Pep","PMA"))
FP$Condition <- as.character("Frozen")

Fresh<-read.csv("Fresh.csv")
splStudyID = strsplit(as.character(Fresh$Sample),":");
Fresh$Sample = as.factor(sapply(splStudyID,"[",2));
splStudyID = strsplit(as.character(Fresh$Sample),"_");
Fresh$Donor = as.factor(sapply(splStudyID,"[",1));
Fresh$Antigen = as.factor(sapply(splStudyID,"[",2));
Fresh$Stim = factor(Fresh$Antigen,levels=c("UN", "WCL","Pep","PMA"))
Fresh$Condition <- as.character("Fresh")

#plot TB
TB<-rbind.fill(Fresh,FP)
TB<-data.table(TB)
TB<-TB[!(Donor=="HD136")]
TB<-TB[!(Donor=="HD398")]
TB<-TB[!(Donor=="HD118")]
TB<-TB[!(Stim=="PMA")]
ggplot(TB, aes(x=Stim, y=CD4.TH1, color=FP)) + 
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of TH1 Cytokine+ CD4+ T cells") +
  ggtitle("TH1") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")
ggplot(TB, aes(x=Stim, y=CD4.TH1.2, color=FP)) + 
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of TH1/2 Cytokine+ CD4+ T cells") +
  ggtitle("TH1/TH2") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")
ggplot(TB, aes(x=Stim, y=CD4.TH2, color=FP)) + 
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of TH2 Cytokine+ CD4+ T cells") +
  ggtitle("TH2") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

#####summary plots
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

pd <- position_dodge(0.1) # move them .05 to the left and right

###
Dat1<-summarySE(FP, measurevar="CD4.TH1", groupvars=c("Stim","FP"))
ggplot(Dat1, aes(x=Stim, y=CD4.TH1, color=FP)) + 
  geom_errorbar(aes(ymin=CD4.TH1-se, ymax=CD4.TH1+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of TH1 Cytokine+ CD4+ T cells") +
  ggtitle("TH1") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

Dat2<-summarySE(FP, measurevar="CD4.TH1.2", groupvars=c("Stim","FP"))
ggplot(Dat2, aes(x=Stim, y=CD4.TH1.2, color=FP)) + 
  geom_errorbar(aes(ymin=CD4.TH1.2-se, ymax=CD4.TH1.2+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of TH1/2 Cytokine+ CD4+ T cells") +
  ggtitle("TH1/TH2") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

Dat3<-summarySE(FP, measurevar="CD4.TH2", groupvars=c("Stim","FP"))
ggplot(Dat3, aes(x=Stim, y=CD4.TH2, color=FP)) + 
  geom_errorbar(aes(ymin=CD4.TH2-se, ymax=CD4.TH2+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of TH2 Cytokine+ CD4+ T cells") +
  ggtitle("TH2") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")

#####Oregon Green
OG<-read.csv("OG.csv")
splStudyID = strsplit(as.character(OG$Sample),":");
OG$Sample = as.factor(sapply(splStudyID,"[",2));
splStudyID = strsplit(as.character(OG$Sample),"_");
OG$Donor = as.factor(sapply(splStudyID,"[",1));
OG$Antigen = as.factor(sapply(splStudyID,"[",2));
OG$Stim = factor(OG$Antigen,levels=c("UN-UN","UN-PMA","WCL-UN", "WCL-PMA","SEB-UN","SEB-PMA"))
OG$Time<-as.character(OG$Time)

Dat4<-summarySE(OG, measurevar="CD4_OG.L", groupvars=c("Stim","Time"))
ggplot(Dat4, aes(x=Stim, y=CD4_OG.L, color=Time)) + 
  geom_errorbar(aes(ymin=CD4_OG.L-se, ymax=CD4_OG.L+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of OG Low CD4+ T cells") +
  ggtitle("Oregon Green Proliferation") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  theme(legend.position="bottom")

Dat5<-summarySE(OG, measurevar="CD4_Cytokine", groupvars=c("Stim","Time"))
ggplot(Dat5, aes(x=Stim, y=CD4_Cytokine, color=Time)) + 
  geom_errorbar(aes(ymin=CD4_Cytokine-se, ymax=CD4_Cytokine+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of Cytokine+ CD4+ T cells") +
  ggtitle("Cytokine Production") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  theme(legend.position="bottom")

Dat6<-summarySE(OG, measurevar="CD8_OG.L", groupvars=c("Stim","Time"))
ggplot(Dat6, aes(x=Stim, y=CD8_OG.L, color=Time)) + 
  geom_errorbar(aes(ymin=CD8_OG.L-se, ymax=CD8_OG.L+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of OG Low CD8+ T cells") +
  ggtitle("Oregon Green Proliferation") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  theme(legend.position="bottom")

Dat7<-summarySE(OG, measurevar="CD8_Cytokine", groupvars=c("Stim","Time"))
ggplot(Dat7, aes(x=Stim, y=CD8_Cytokine, color=Time)) + 
  geom_errorbar(aes(ymin=CD8_Cytokine-se, ymax=CD8_Cytokine+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of Cytokine+ CD8+ T cells") +
  ggtitle("Cytokine Production") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  theme(legend.position="bottom")

#stats
melted<-melt(FP, id.vars=c("Stim", "FP"))
melted$value<-as.numeric(as.character(melted$value))
stats<-as.data.frame(
  ddply(melted, c("Stim", "FP", "variable"), summarise,
        mean = mean(value), sd = sd(value),
        sem = sd(value)/sqrt(length(value))))


###OG Antigen Stims
OG<-read.csv("OG Antigen Stim.csv")
splStudyID = strsplit(as.character(OG$Sample),":");
OG$Sample = as.factor(sapply(splStudyID,"[",2));
splStudyID = strsplit(as.character(OG$Sample),"_");
OG$Donor = as.factor(sapply(splStudyID,"[",1));
OG$Antigen = as.factor(sapply(splStudyID,"[",2));
splStudyID = strsplit(as.character(OG$Antigen)," ");
OG$Antigen = as.factor(sapply(splStudyID,"[",1));
OG$Stim = factor(OG$Antigen,levels=c("UN","WCL","PPD", "SEA","SWAP","SEB"))


Dat8<-summarySE(OG, measurevar="CD4.OG.L", groupvars="Stim")
ggplot(Dat8, aes(x=Stim, y=CD4.OG.L)) + 
  geom_errorbar(aes(ymin=CD4.OG.L-se, ymax=CD4.OG.L+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=10) +
  xlab("Stimulation") +
  ylab("Frequency of OG Low CD4+ T cells") +
  ggtitle("Oregon Green Proliferation") +
  theme_bw()+
  theme(text = element_text(size = 30))+
  theme(legend.position="bottom")