# All of this generates a fake data table that contains a "key" variable (i.e. donor)
# a condition variable (i.e. "stim") and a set of columns of output data
# the goal is for every donor, in each column of output data subtract out the value in that column
# that corresponds to a specific value in the condition column
donor<-c(rep("john", 6), rep("mary", 6), rep("jack", 6), rep("suzy", 6))
stim<-rep(c("un","w","p","s","sw","sb"),4)
ifng<-rnorm(24, 2)
tnfa<-rnorm(24,10)
IL4<-rnorm(24,5)
IL5<-rnorm(24,7)
IL10<-rnorm(24,13)
IL13<-rnorm(24,19)
IL17<-rnorm(24,17)
IL21<-rnorm(24,11)
IL22<-rnorm(24,3)

fake<-data.frame(cbind(donor,stim,ifng,tnfa,IL4, IL5, IL10, IL13, IL17, IL21, IL22))
remove(donor,stim,ifng,tnfa,IL4, IL5, IL10, IL13, IL17, IL21, IL22)
#can only use data frames to use bracketed string column indexing so make sure this is a data frame not a data table

#when binding data together the output data becomes facgtors and we need them to be numeric
class_change<-function(datatable){
    y<-dplyr::select(datatable, donor, stim)
    x<-dplyr::select(datatable, -donor, -stim)
    for(i in 1:length(x)){x[[i]]<-as.numeric(as.character(x[[i]]))}
    df<-cbind(y, x)
}
fake<-class_change(fake)
#we want to copy the original data table to reference back to and check our math
fakecopy<-fake

#This is going to generate the list of columns names that we want to iterate over
numeric.only <- function(X,...){
    returnCols <- names(X)
    a<-sapply(X, is.numeric)
    print(returnCols[a == "TRUE"])
}
colnms<-numeric.only(fake)

dt<-split(fake, fake[,"donor"]) #generate a list of donor specific mini data frames 

#split --> lists in alphabetical order not by original key order
#To solve this I arranged the original data frame in alphabetical order

fake<-arrange(fake, donor)
#can we instead arrange the split to match the original order

#for each column name in the list of column names
for(col in colnms){
    #this is a check to make sure it goes through all the columns
    print(col) 
    #dt.r will make a list of length of dt where each item in the list is the new column for a particular donor
    dt.r<-lapply(dt, function(g)
        #take the values in each column and subtract out the value indicated by an "un" the stim column 
        (g[,col] - subset(g, g$stim=="un")[,col]))
    #unlist the list so you get a single column and insert it back into the original data frame
    fake[,col]<-unlist(dt.r)    
}


