PBMC<-read.csv("PBMC.csv")
PBMC$Viability<-as.numeric(PBMC$Viability)
PBMC$SM<-factor(PBMC$SM, levels=c("Negative","Schisto+","HEL"))

ggplot(PBMC, aes(x=TB.Status, y=Viability, color=SM))+
  xlab("") +
  ylab("Frequency of Viable Cells") +
  ggtitle("Viability") +
  geom_boxplot(outlier.shape = NA, size=1.8, position=dodge)+
  scale_color_manual(values=Infection, name="Helminth Status")+
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=24,face="bold"),plot.title=element_text(size=30,face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.key.size = unit(2.5, "cm"))+
  theme(legend.text = element_text(size = 20))

ggplot(PBMC, aes(x=TB.Status, y=Cells.Vial, color=SM))+
  xlab("") +
  ylab("Number Cells Per Vial Thawed") +
  ggtitle("Cell Counts") +
  geom_boxplot(outlier.shape = NA, size=1.8, position=dodge)+
  scale_color_manual(values=Infection, name="Helminth Status")+
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=24,face="bold"),plot.title=element_text(size=30,face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.key.size = unit(2.5, "cm"))+
  theme(legend.text = element_text(size = 20))+
  coord_cartesian(ylim = c(0,2E7))

#Boolean plots
ggplot(DF, aes(x=label, y=value, color=SM))+
  ylab("Frequency Cytokine+ CD4+ T-cells (%)") +
  ggtitle(DF$TB) +
  geom_boxplot(outlier.shape = NA, size=1.8, position=dodge)+
  scale_color_manual("Helminth Status", values=Infection)+
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=24,face="bold"),plot.title=element_text(size=30,face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.key.size = unit(2.5, "cm"))+
  theme(legend.text = element_text(size = 20))+
  coord_cartesian(ylim=c(0,25))

#reg cytokine

scale_color_manual(values=c("#009999","#663399","#CC0066"), name="Helminth Status")+




