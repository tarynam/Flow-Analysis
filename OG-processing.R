library(dplyr)
library(data.table)

setwd("/Applications/Old Computer/Day Lab/Flow-Data/oregon-green")
filenames = list.files(pattern="csv")
for(i in 1:length(filenames)){
    print(filenames[i])
    assign(filenames[i], read.table(filenames[i], header=TRUE, sep=";"))
}

l.df<-lapply(ls(pattern=".csv"), function(x) get(x))
combined<-dplyr::bind_rows(l.df)
fwrite(combined, file="/Applications/Old Computer/Day Lab/Flow-Data/Master-OG2.csv")

