library(dplyr)

setwd("/Users/tarynam/Desktop/ICS-CD4")
filenames = list.files(pattern="csv")
for(i in 1:length(filenames)){
    print(filenames[i])
    assign(filenames[i], read.table(filenames[i], header=TRUE, sep=";"))
}

l.df<-lapply(ls(pattern=".csv"), function(x) get(x))
combined<-dplyr::bind_rows(l.df)

#this has to be edited still 
combined$Sample<-paste(combined$donor, combined$Stim, sep="_")
remove<-c("NK2126_HC_X_SEA", 
          "NK2450_HC_SM_SEA" , 
          "NK2402_LTBI_SM_SWAP",
          "NK2342_TB_X_SWAP",
          "NK2136_TB_X_WCL",
          "NK2168_HC_SM_WCL",
          "NK2421_LTBI_SM_PEP",
          "NK2063_TB_X_SEA",
          "NK2115_LTBI_X_PEP",
          "NK2115_LTBI_X_WCL",
          "NK2171_HC_X_WCL",
          "NK2162_HC_X_SEA")
combined<-dplyr::filter(combined, !(Sample %in% remove))
write.csv(combined, "/Users/tarynam/Desktop/ICS-CD4/clean-data/Master-ICS.csv")


