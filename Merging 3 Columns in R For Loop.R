TB.1<-c("LTBI", "LTBI", "LTBI", "LTBI", "LTBI", 
        "ACTIVE", "ACTIVE", "ACTIVE", "HEALTHY CONTROL", "HEALTHY CONTROL",
        NA,NA,NA,NA,NA)
TB.2<-c(NA, NA, "LTBI", "LTBI", "HEALTHY CONTROL", 
        "ACTIVE", "ACTIVE", NA, "HEALTHY CONTROL", "HEALTHY CONTROL",
        "HEALTHY CONTROL","ACTIVE","LTBI",NA,NA)
TB.3<-c(NA, "ACTIVE", NA, NA, "HEALTHY CONTROL", 
        NA, NA, NA,"HEALTHY CONTROL", NA,
        "HEALTHY CONTROL","ACTIVE","LTBI","LTBI", "HEALTHY CONTROL")
library(data.table)
TB<-as.data.table(cbind(TB.1,TB.2,TB.3))
TB[,"TB.Final"]<-NA

###I'm trying to write a for-loop that combines the information in columns 1-3
#such that if TB.1 is not NA then TB.Final takes that value
#if TB.1 is NA then TB.Final takes the value in TB.2
#if TB.2 is NA then TB.Final takes the value in TB.3

#This one doesn't incorporate the information in TB.3 without an error message
for (i in 1:length(TB$TB.Final)){
        if(!is.na(TB$TB.1[i])){
                TB$TB.Final[i] <- TB$TB.1[i]}
        
        else if(!is.na(TB$TB.2[i])){
                TB$TB.Final[i] <- TB$TB.2[i]} 
        
        else{
                TB$TB.Final[i] <- TB$TB.3[i]} 
}

#This one also doesn't incorporate the information in TB.3 without an error message
for (i in 1:length(TB$TB.Final)){
        if(!is.na(TB$TB.1[i])){
                TB$TB.Final[i] <- TB$TB.1[i]}
        
        else if(is.na(TB$TB.1[i])){
                TB$TB.Final[i] <- TB$TB.2[i]} 
        
        else if (is.na(TB$TB.1[i] & is.na(TB$TB.2[i]))){
                TB$TB.Final[i] <- TB$TB.3[i]} 
}


