setwd("/Applications/Old Computer/Day Lab/Flow-Data/clean")
filenames =list.files(pattern="OG")

pbmc<-read.csv("DonorInfo_clean.csv")
pbmc<-dplyr::select(pbmc, Donor, TB, SM, used, OG_Excluded)%>%
    dplyr::filter(grepl("OG", pbmc$used))

test<-read.csv("OG_GDonly_OG_summary_clean.csv")
#Selecting and renaming variables
test2<-dplyr::select(test, Donor, Stim, TB, SM)

inter<-intersect(pbmc$Donor, test$Donor)
pbmc$Donor[!(pbmc$Donor %in% inter)]
test$Donor[!(test$Donor %in% inter)]

compare<-merge(test2, pbmc, by="Donor", all=TRUE)

for(i in 1:length(compare)){compare[[i]]<-as.character(compare[[i]])}
length(which(compare$TB.x != compare$TB.y))
length(which(compare$SM.x != compare$SM.y))



