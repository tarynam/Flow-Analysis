names<-c("ID1","ID2","ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10")
data<-matrix(data=rep(1:8, 5), nrow=10)
data<-as.data.table(cbind(names,data))
names(data)<-c("ID","CXCR3_MFI","blank","blank2","CXCR3_Freq")

#I'm super lost on how to combine various tools to achieve the following
#I want to make certain values of columns that contain the string "CXCR3" NA
#based on the ID in the ID column

#I know I can select columns containing "CXCR3" using contains in dplyr
library(dplyr)
df1<-data %>% select(contains("CXCR3"))
#but I can't combine that with "ID" these both give errors
df1<-data %>% select(contains("CXCR3|ID"))
df1<-data %>% select(contains("CXCR3") | contains("ID"))

#So then I tried a million variations of grep commands
#This one gives me a 3 rows all columns
df2 <- data [grep("CXCR3|ID", names(data))]
#This one has an error
df3 <- data [grepl("CXCR3|ID", names(data))]
#This one give me a logical vector of all FALSE
df4 <- grepl("ID | CXCR3", names(data))

#this gives me the correct column names but I can never seem to subset the data based on this
colnames<-names(data)[grep("ID|CXCR3", names(data))]
df<-select(data, (names(data) %in% colnames))
df<-data[(names(data) %in% colnames)]
df<-subset(data, names(data) %in% colnames)

#So then I gave up but this is what I want to do 
#but without having to type out every column name <-NA
important_cols = grep('CXCR3', names(data))
for (i in 1:length(data)){
        if(data$ID[i] %in% c("ID1", "ID4", "ID5")){
                for(j in important_cols){
                        data[i,j] <- NA;
                }        
        }
}


