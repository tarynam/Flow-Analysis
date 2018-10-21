DF<-read.csv("/Applications/Old Computer/Day Lab/Data_Flow/Counts.csv", stringsAsFactors = TRUE)
library(dplyr)
DF<-filter(DF,!(Sample=="Mean")) %>%
  filter(!Sample=="StdDev") %>%
  filter(!Sample=="Sample") %>%
  filter(!Sample=="NK2126_HC_X_SEA_005.fcs") %>%
  filter(!Sample=="NK2450_HC_SM_SEA_047.fcs") %>%
  filter(!Sample=="NK2402_LTBI_SM_SWAP_024.fcs") %>%
  filter(!Sample=="NK2342_TB_X_SWAP_006.fcs") %>%
  filter(!Sample=="NK2136_TB_X_WCL_003.fcs") %>%
  filter(!Sample=="NK2168_HC_SM_WCL_040.fcs") %>%
  filter(!Sample=="NK2421_LTBI_SM_PEP_010.fcs") %>%
  filter(!Sample=="NK2063_TB_X_SEA_005.fcs") %>%
  filter(!Sample=="NK2115_LTBI_X_PEP_010.fcs") %>%
  filter(!Sample=="NK2115_LTBI_X_WCL_009.fcs") %>%
  filter(!Sample=="NK2171_HC_X_WCL_009.fcs") %>%
  filter(!Sample=="NK2162_HC_X_SEA_023.fcs")

source("/Applications/Old Computer/DataScience/Flow Shit/tryna function.R")
library(data.table)
data<-Split(DF, column="Sample", splitter="_", one="ID", two="TB", three="SM", four="Stim")

DF2<-read.csv("/Applications/Old Computer/Day Lab/Data_Flow/Clean/Cytokine.csv", stringsAsFactors = TRUE)
DF2<-DF2[,c(3,24:27)]
df1<-merge(data, DF2, by="Sample")

df2<-df1[,2:17]
colnames(df2)<-compass_nms
data2<-cbind(df1[,1], df1[,18:21], df2, df1[,22:25])

write.csv(data2,"/Applications/Old Computer/Day Lab/Data_Flow/Clean/All Counts.csv")




