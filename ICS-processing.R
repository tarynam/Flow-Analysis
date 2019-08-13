library(dplyr)

setwd("/Applications/Old Computer/Day Lab/Flow-Data/ICS-CD4")
filenames = list.files(pattern="csv")
for(i in 1:length(filenames)){
    print(filenames[i])
    assign(filenames[i], read.table(filenames[i], header=TRUE, sep=";"))
}

NK2054_TB_SM.csv$donor<-"NK2054_TB_SM"

l.df<-lapply(ls(pattern=".csv"), function(x) get(x))
combined<-dplyr::bind_rows(l.df)
write.csv(combined, "/Applications/Old Computer/Day Lab/Flow-Data/Master-ICS.csv")

