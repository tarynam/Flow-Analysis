library(dplyr)

setwd("/Users/tarynam/Desktop/ICS-CD4")
filenames = list.files(pattern="csv")
for(i in 1:length(filenames)){
    print(filenames[i])
    assign(filenames[i], read.table(filenames[i], header=TRUE, sep=";"))
}

l.df<-lapply(ls(pattern=".csv"), function(x) get(x))
combined<-dplyr::bind_rows(l.df)
write.csv(combined, "/Users/tarynam/Desktop/ICS-CD4/clean-data/Master-ICS.csv")


