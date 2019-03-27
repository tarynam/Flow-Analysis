#used to generate a list of character strings that are the names of the columns containing numeric data
numeric.only <- function(X,...){
    returnCols <- names(X)
    a<-sapply(X, is.numeric)
    print(returnCols[a == "TRUE"])
}

#takes a dataframe (v. important its a dataframe), splitvar is a character string of a column name
#colnms is a list of strings indicating which columns to subtract background from
#condition is a character string with the key to indicate what is being subtracted
background.subtract<-function(df, splitvar, keyvar, condition){
    dt<-split(df, df[, splitvar]) #generate a list of mini data frames for each level of splitvar
    #split --> lists in alphabetical order not by original key order
    #To solve this I arranged the original data frame in alphabetical order
    #can we instead arrange the split to match the original order
    
    colnms<-numeric.only(df)
    
    #for each column name in the list of column names
    for(col in colnms){
        #this is a check to make sure it goes through all the columns
        print(col) 
        #dt.r will make a list of length of dt where each item in the list is the new column for a particular donor
        dt.r<-lapply(dt, function(g)
            #take the values in each column and subtract out the value indicated by an "un" the stim column 
            (g[,col] - subset(g, g[,keyvar]==condition)[,col]))
        #unlist the list so you get a single column and insert it back into the original data frame
        df[,col]<-unlist(dt.r)
    }
    df
}

filenames = c("CD4_Cytokine.csv", "CD8_Cytokine.csv", "GD_Cytokine.csv")
l<-list(CD4_Cytokine.csv, CD8_Cytokine.csv, GD_Cytokine.csv)
l<-lapply(l, function(x) arrange(x, Donor))
l<-lapply(l, function(x) background.subtract(x, "Donor", "Stim", "UN"))
for (i in 1:length(filenames)) {
    assign(filenames[i], 
           data.frame(l[i]))}
CD4_Cytokine.csv[CD4_Cytokine.csv<0]<-0
CD8_Cytokine.csv[CD8_Cytokine.csv<0]<-0
GD_Cytokine.csv[GD_Cytokine.csv<0]<-0