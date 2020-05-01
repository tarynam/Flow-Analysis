library(dplyr)
library(data.table)

setwd("/Users/tarynam/Desktop/Lisa/")
compass<-read.csv("MIMI_compass.csv")
tiger<-dplyr::select(compass, Sample, unique, ID, Stim) #this should just me the meta data
kitten<- compass[,2:17] #This should just be the count columns

### TWO WAYS TO RELABEL CYTOKINES 

# IF YOU HAVE USE THE .P .N NOMENCLATURE
# do this for each cytokine make sure the capitalization matches whatever you have
names(kitten) <- gsub("IL2.p", "&IL2", names(kitten))
names(kitten) <- gsub("IL2.n", "&!IL2", names(kitten))
# Whatever the first cytokine is gets just a "!" instead of a "&!"
names(kitten) <- gsub("&!IFNG", "!IFNG", names(kitten))

# IF TAKING THE NAMES STRAIGHT FROM FLOWJO AND NOT MANUALLY EDITING THE FILE
# These are the original column names in order. 
# The nice thing is your can manually reorder the cytokines for the future graphs
# Ex: I flipped the TNF and IL17 so that the cytokine ordering looks better in the plots
#i+2+17+t+  --->  IFNγ&IL2&TNFα&IL17
#i+2+17+t-  --->  IFNγ&IL2&!TNFα&IL17
#i+2+17-t+  --->  IFNγ&IL2&TNFα&!IL17
#i+2+17-t-  --->  IFNγ&IL2&!TNFα&!IL17
#i+2-17+t+  --->  IFNγ&!IL2&TNFα&IL17
#i+2-17+t-  --->  IFNγ&!IL2&!TNFα&IL17
#i+2-17-t+  --->  IFNγ&!IL2&TNFα&!IL17
#i+2-17-t-  --->  IFNγ&!IL2&!TNFα&!IL17
#i-2+17+t+  --->  !IFNγ&IL2&TNFα&IL17
#i-2+17+t-  --->  !IFNγ&IL2&!TNFα&IL17
#i-2+17-t+  --->  !IFNγ&IL2&TNFα&!IL17
#i-2+17-t-  --->  !IFNγ&IL2&!TNFα&!IL17
#i-2-17+t+  --->  !IFNγ&!IL2&TNFα&IL17
#i-2-17+t-  --->  !IFNγ&!IL2&!TNFα&IL17
#i-2-17-t+  --->  !IFNγ&!IL2&TNFα&!IL17
#i-2-17-t-  --->  !IFNγ&!IL2&!TNFα&!IL17

newnames<-c("IFNγ&IL2&TNFα&IL17",
             "IFNγ&IL2&!TNFα&IL17",
             "IFNγ&IL2&TNFα&!IL17",
             "IFNγ&IL2&!TNFα&!IL17",
             "IFNγ&!IL2&TNFα&IL17",
             "IFNγ&!IL2&!TNFα&IL17",
             "IFNγ&!IL2&TNFα&!IL17",
             "IFNγ&!IL2&!TNFα&!IL17",
             "!IFNγ&IL2&TNFα&IL17",
             "!IFNγ&IL2&!TNFα&IL17",
             "!IFNγ&IL2&TNFα&!IL17",
             "!IFNγ&IL2&!TNFα&!IL17",
             "!IFNγ&!IL2&TNFα&IL17",
             "!IFNγ&!IL2&!TNFα&IL17",
             "!IFNγ&!IL2&TNFα&!IL17",
             "!IFNγ&!IL2&!TNFα&!IL17")    
names(kitten)<-newnames

## PHA
#The donor names in the unstim data and the stim data have to match
#this doesn't always happen when we remove specific samples
#so this is a way to make sure that they match by getting the "intersecting" donors
unstim_names<-subset(tiger$unique,tiger$Stim=="unstim")
stim_names<-subset(tiger$unique,tiger$Stim=="PHA")
inter<-intersect(unstim_names , stim_names)

#Now we use that list of donor IDs to create a subsetted count dataframe for the unstim
lion<-subset(kitten, tiger$Stim =="unstim" & tiger$unique %in% inter)
#Need to give the data table rownames so we can match unstim to stim
rownames(lion) <- tiger$unique[(tiger$Stim =="unstim" & tiger$unique %in% inter)]
#repeat for the stim
bobcat<-subset(kitten, tiger$Stim =="PHA" & tiger$unique %in% inter)
rownames(bobcat) <- tiger$unique[(tiger$Stim =="PHA" & tiger$unique %in% inter)]
tiger_s<-data.frame(subset(tiger,tiger$Stim=="PHA"))

fit_PHA<-COMPASS::SimpleCOMPASS (n_s=bobcat, #the unstim count data frame
                                 n_u=lion, #the stim count data frame
                                 meta=tiger_s, #the metadata
                                 individual_id = "unique", #the way to match the meta data to the count data
                                 iterations = 10000, replications = 8, #math stuff
                                 verbose = TRUE) #print info as it runs
#Sometimes the rownames aren't preserved
#you may need to do use the next line which I have commented out
#rownames(fit_PHA$fit$gamma)<-rownames(bobcat)

#save the fit and the scores for the future
save(fit_PHA, file="fit_PHA.RData")
write.csv(scores(fit_PHA), 
          "compass_MIMI_PHA.csv")

##you can plot the compass data immediately or later all at once
png("name your plot here.png")
plot(fit_PHA, 
     order_by_max_functionality = FALSE, #Doesn't rearrange columns by scores
     markers=c("IL17","TNFα","IL2","IFNγ"), #specify the order of the rows IFNg is the top
     row_annotation=c("TB","SM"), #the column names in tiger to annotate rows
     border_color=NA,
     fontsize=14, 
     threshold=-1) #Doesn't remove columns with low scores 
#There are a lot more options for specifying things but we'll cross that bridge later
dev.off()





