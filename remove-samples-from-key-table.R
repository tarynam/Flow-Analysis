# All of this generates a fake data table that contains a "key" variable (i.e. donor)
# a condition variable (i.e. "stim") and a set of columns of output data
# the goal is for every donor, in each column of output data subtract out the value in that column
# that corresponds to a specific value in the condition column
donor<-c(rep("john", 6), rep("mary", 6), rep("jack", 6), rep("suzy", 6))
stim<-rep(c("un","w","p","s","sw","sb"),4)
Cyt<-sample(20:200, 24)
TH1<-sample(20:200, 24)
TH2<-sample(20:200, 24)
TH1.2<-sample(20:200, 24)

Cyt_ifng<-sample(100:1000, 24)
Cyt_tnfa<-sample(100:1000, 24)
Cyt_IL4<-sample(100:1000, 24)
Cyt_IL13<-sample(100:1000, 24)
TH1_ifng<-sample(100:1000, 24)
TH1_tnfa<-sample(100:1000, 24)
TH1_IL4<-sample(100:1000, 24)
TH1_IL13<-sample(100:1000, 24)
TH2_ifng<-sample(100:1000, 24)
TH2_tnfa<-sample(100:1000, 24)
TH2_IL4<-sample(100:1000, 24)
TH2_IL13<-sample(100:1000, 24)
TH1.2_ifng<-sample(100:1000, 24)
TH1.2_tnfa<-sample(100:1000, 24)
TH1.2_IL4<-sample(100:1000, 24)
TH1.2_IL13<-sample(100:1000, 24)

fake_count<-data.frame(cbind(donor,stim, Cyt, TH1, TH2, TH1.2))
fake_mfi<-data.frame(cbind(donor,stim,Cyt_ifng,Cyt_tnfa,Cyt_IL4,Cyt_IL13,TH1_ifng,TH1_tnfa,
                           TH1_IL4,TH1_IL13,TH2_ifng,TH2_tnfa,TH2_IL4,TH2_IL13,TH1.2_ifng,TH1.2_tnfa,
                           TH1.2_IL4,TH1.2_IL13))
remove(donor, stim, Cyt, TH1, TH2, TH1.2,Cyt_ifng,Cyt_tnfa,Cyt_IL4,Cyt_IL13,TH1_ifng,TH1_tnfa,
       TH1_IL4,TH1_IL13,TH2_ifng,TH2_tnfa,TH2_IL4,TH2_IL13,TH1.2_ifng,TH1.2_tnfa,
       TH1.2_IL4,TH1.2_IL13)
#can only use data frames to use bracketed string column indexing so make sure this is a data frame not a data table

#when binding data together the output data becomes facgtors and we need them to be numeric
class_change<-function(datatable){
    y<-dplyr::select(datatable, donor, stim)
    x<-dplyr::select(datatable, -donor, -stim)
    for(i in 1:length(x)){x[[i]]<-as.numeric(as.character(x[[i]]))}
    df<-cbind(y, x)
}
fake_count<-class_change(fake_count)
fake_mfi<-class_change(fake_mfi)


#brute force reshape strategy
melt<-melt(fake_count, id.vars=c("donor", "stim"))
melt<-dplyr::rename(melt, count=value, cell_subset=variable)
melt_mfi<-melt(fake_mfi, id.vars=c("donor","stim"))
splitter<-function(data,column,splitter,one=NULL,two=NULL){
    library(data.table)
    
    labels<-tstrsplit(data[,column], splitter)
    data[,one]<-labels[[1]]
    data[,two]<-labels[[2]]
    return(data)
}
melt_mfi<-splitter(melt_mfi, "variable", "_", one="cell_subset", two="attribute")
melt_mfi<-select(melt_mfi, donor, stim, cell_subset, attribute, value)
test<-full_join(melt, melt_mfi)
for(r in 1:nrow(test)){
    if(test$count[r]<50){
        test$value[r]<-NA
    }
}
test2<-spread(test, attribute,value)

###try with a double for loop
for(row in 1:nrow(fake_mfi)){
    keep <- rownames(fake_count)[fake_count$Cyt<50]
    if(row %in% keep){
        fake_mfi$Cyt_ifng[row]<-NA
    }
    }
