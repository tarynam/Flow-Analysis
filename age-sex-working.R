library(ggplot2)
library(ggpubr)
library(dplyr)

DonorInfo<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/clean/ICS_Donor_Info_clean.csv")
DonorInfo<-filter(DonorInfo, SM %in% c("X", "SM"))
DonorInfo$SM<-factor(DonorInfo$SM, levels=c("X", "SM"))
DonorInfo$TB<-factor(DonorInfo$TB, levels=c("HC", "LTBI", "TB"))
DonorInfo$disease<-paste(DonorInfo$TB, DonorInfo$SM, sep="_")
DonorInfo$disease<-factor(DonorInfo$disease, levels=c("HC_X", "HC_SM",
                                                      "LTBI_X", "LTBI_SM",
                                                      "TB_X", "TB_SM"))
my_comparisons <- list(c("HC_SM", "HC_X"), c("LTBI_SM", "LTBI_X"), c("TB_SM", "TB_X"), 
                       c("HC_SM", "LTBI_SM"), c("TB_SM", "LTBI_SM"), c("HC_SM", "TB_SM"),
                       c("HC_X", "LTBI_X"), c("TB_X", "LTBI_X"), c("HC_X", "TB_X"))

IQR_vals <- function(datatable, column) {
    library(data.table)
    y<-split(datatable, datatable$disease)
    x<-lapply(y, function(y) quantile(y[,column], probs=c(0.5, 0.25, 0.75), na.rm=TRUE))
    #quantile is a stats function that takes a vector "x" and returns the value at a specified quantile
    data.frame(t(x))
}
IQR_pvals<-function(data, column){
    x<-split(data,data$TB)
    #does a wilcox test between SM- and SM+ for whatever column you call in the column argument
    A<-lapply(x, function(g) wilcox.test(g[,column]~g[,"SM"], exact=FALSE, na.action = na.pass))
    pvals<-c(A$HC$p.value, A$LTBI$p.value, A$TB$p.value)
    fdr_adj_pvals<-round(p.adjust(pvals, method="fdr"), 4)
    t(fdr_adj_pvals)         
}


###N
n<-t(DonorInfo %>% dplyr::count(disease))
names(n)<-n[1,]

## AGE
fit<-lm(age~disease, data=DonorInfo)
total.p<-anova(fit)$`Pr(>F)`
#Age is different
ggplot(DonorInfo, aes(x=disease, y=age)) +
    geom_boxplot(aes(fill=TB, alpha=SM), size=1, position=position_dodge(width = 1), outlier.shape=NA) +
    geom_jitter(width=.1,height=0, shape=16,size=4)+
    scale_alpha_manual(values=c(0.5,1))+
    scale_fill_manual(values = c("#1a9850" , "#2166ac", "#b2182b"))+
    theme_classic()+ theme(legend.position="none")+
    theme(text = element_text(size=20), axis.text.x = element_text(angle=90, vjust=0.6)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(comparisons = my_comparisons,  p.adjust="bonferroni", aes(label=..p.adj..))+
    labs(subtitle = "These aren't corrected for whatever reason will figure out later")
AGE<-cbind(IQR_vals(DonorInfo, "age"),
           cbind(IQR_pvals(DonorInfo, "age"), total.p))
AGE<-dplyr::rename(AGE, HC.p = V1, LTBI.p = V2, TB.p = V3)
AGE$parameter<-"Age"

## SEX
female<-(DonorInfo %>% dplyr::group_by(sex) %>% dplyr::count(disease))[1:6,2:3]
male<-(DonorInfo %>% dplyr::group_by(sex) %>% dplyr::count(disease))[7:12,2:3]
df.sex<-merge(female, male, by="disease")
df.sex<-dplyr::rename(df.sex, female=n.x, male=n.y)
df.sex$femfreq = paste("(",
    round(df.sex$female/(df.sex$female+df.sex$male),2)*100, 
    "%", ")", sep="")
df.sex$malefreq = paste("(",
    round(df.sex$male/(df.sex$female+df.sex$male),2)*100,
    "%", ")", sep="")

total.p<-round(chisq.test(df.sex[,2:3])$p.value, 3)
HC.p<-round(chisq.test(df.sex[1:2,2:3])$p.value, 3)
LTBI.p<-round(chisq.test(df.sex[3:4,2:3])$p.value, 3)
TB.p<-round(chisq.test(df.sex[5:6,2:3])$p.value, 3)
pvals<-data.frame(cbind(HC.p, LTBI.p, TB.p, total.p))
female<-paste(df.sex$female, df.sex$femfreq, sep=" ")
male<-paste(df.sex$male, df.sex$malefreq, sep=" ")
d1<-data.frame(rbind(female, male))
names(d1)<-df.sex$disease
SEX<-cbind(d1, pvals)
SEX$parameter<-row.names(SEX)

##RACE
RACE<-paste(n[2,], "(100%)", sep=" ")

## QFT
qft<-filter(DonorInfo, TB!="TB")
qft$TB<-factor(qft$TB, levels=c("HC", "LTBI"))
qft$disease<-factor(qft$disease, levels=c("HC_X","HC_SM",
                                          "LTBI_X", "LTBI_SM"))
fit<-lm(QFT~disease, qft)
total.p<-anova(fit)$`Pr(>F)`
QFT<-IQR_vals(qft, "QFT")
x<-split(qft,qft$TB)
#does a wilcox test between SM- and SM+ for whatever column you call in the column argument
A<-lapply(x, function(g) wilcox.test(g$QFT~g$SM, exact=FALSE, na.action = na.pass))
pvals<-c(A$HC$p.value, A$LTBI$p.value, total.p)

## Schisto
eggs<-filter(DonorInfo, SM=="SM")
fit<-lm(SchistosomaIntensity~TB, eggs)
anova(fit)$`Pr(>F)`


####NEED TO DEAL WITH TOTALS


#plots to check age vs
compass<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/clean/ICS_Compass_scores_complete.csv")
CD4<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/clean/CD4_clean.csv")
CD4_cyt<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/clean/CD4_Cytokine_clean.csv")

test<-merge(CD4, DonorInfo, by="Donor")
test<-filter(test, Stim=="UN")
test.age.plot<-function(datatable, column){
    g<-ggplot(datatable, aes(x=age, y=datatable[,column], color=sex))
    g<- g + geom_point()+theme_bw()+geom_smooth(method="lm", show.legend=TRUE)
    g<- g + labs(title=column)
    g
}
listnames<-names(test)
for(l in listnames){
    print(l)
    p<-test.age.plot(test, l)
    print(p)
}
### By Cytokine
test<-merge(CD4_cyt, DonorInfo, by="Donor")
test<-filter(test, Stim=="WCL")%>%
    select(-disease, -Excluded, -Cells.Vial, -Viability, -SM.y,- TB.y, -X.y, -SM.x,- TB.x, -X.x,-Sample)
test.age.plot<-function(datatable, column){
    g<-ggplot(datatable, aes(x=age, y=datatable[,column], color=sex))
    g<- g + geom_point()+theme_bw()+geom_smooth(method="lm", show.legend=TRUE)
    g<- g + labs(title=column, subtitle="WCL")
    g
}
listnames<-names(test)
for(l in listnames){
    print(l)
    p<-test.age.plot(test, l)
    print(p)
}

test<-merge(CD4_cyt, DonorInfo, by="Donor")
test<-filter(test, Stim=="PEP")%>%
    select(-disease, -Excluded, -Cells.Vial, -Viability, -SM.y,- TB.y, -X.y, -SM.x,- TB.x, -X.x,-Sample)
test.age.plot<-function(datatable, column){
    g<-ggplot(datatable, aes(x=age, y=datatable[,column], color=sex))
    g<- g + geom_point()+theme_bw()+geom_smooth(method="lm", show.legend=TRUE)
    g<- g + labs(title=column, subtitle="PEP")
    g
}
listnames<-names(test)
for(l in listnames){
    print(l)
    p<-test.age.plot(test, l)
    print(p)
}

## COMPASS
test<-merge(compass, DonorInfo, by="Donor")
test<-filter(test, Stim=="WCL", celltype=="CD4")
test.age.plot<-function(datatable, column){
    g<-ggplot(datatable, aes(x=age, y=datatable[,column], color=sex))
    g<- g + geom_point()+theme_bw()+geom_smooth(method="lm", show.legend=TRUE)
    g<- g + labs(title=column, subtitle="WCL")
    g
}
test.age.plot(test, "FS")
test.age.plot(test, "PFS")

test<-merge(CD4_cyt, DonorInfo, by="Donor")
test<-filter(test, Stim=="PEP", celltype=="CD4")
test.age.plot<-function(datatable, column){
    g<-ggplot(datatable, aes(x=age, y=datatable[,column], color=sex))
    g<- g + geom_point()+theme_bw()+geom_smooth(method="lm", show.legend=TRUE)
    g<- g + labs(title=column, subtitle="PEP")
    g
}
listnames<-names(test)
for(l in listnames){
    print(l)
    p<-test.age.plot(test, l)
    print(p)
}

test.age.plot(test, "FS")
test.age.plot(test, "PFS")


