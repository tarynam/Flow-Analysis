#source("http://www.bioconductor.org/biocLite.R")
#biocLite("COMPASS")
library(COMPASS)
library(dplyr)
library(data.table)
CD4_compass<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_CD4_Counts_clean.csv")
CD4_compass$celltype<-"CD4"
CD8_compass<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_CD8_Counts_clean.csv")
CD8_compass$celltype<-"CD8"
GD_compass<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_GD_Counts_clean.csv")
GD_compass$celltype<-"GD"


test<-merge(CD4_compass, CD8_compass, 
            by=c("Sample", "Donor", "TB", "SM", "Stim"), all=TRUE)

test2<-merge(test, GD_compass, 
             by=c("Sample", "Donor", "TB", "SM", "Stim"), all=TRUE)

CD3<-dplyr::select(test2, Sample, Donor, TB, SM, Stim)

CD3$G_4_13_T<-test2$G_4_13_T + test2$G_4_13_T.x + test2$G_4_13_T.y
CD3$G_4_13_x<-test2$G_4_13_x + test2$G_4_13_x.x + test2$G_4_13_x.y
CD3$G_4_x_T<-test2$G_4_x_T   + test2$G_4_x_T.x  + test2$G_4_x_T.y
CD3$G_4_x_x<-test2$G_4_x_x   + test2$G_4_x_x.x  + test2$G_4_x_x.y

CD3$G_x_13_T<-test2$G_x_13_T + test2$G_x_13_T.x + test2$G_x_13_T.y
CD3$G_x_13_x<-test2$G_x_13_x + test2$G_x_13_x.x + test2$G_x_13_x.y
CD3$G_x_x_T<-test2$G_x_x_T   + test2$G_x_x_T.x  + test2$G_x_x_T.y
CD3$G_x_x_x<-test2$G_x_x_x   + test2$G_x_x_x.x  + test2$G_x_x_x.y

CD3$x_4_13_T<-test2$x_4_13_T + test2$x_4_13_T.x + test2$x_4_13_T.y
CD3$x_4_13_x<-test2$x_4_13_x + test2$x_4_13_x.x + test2$x_4_13_x.y
CD3$x_4_x_T<-test2$x_4_x_T   + test2$x_4_x_T.x  + test2$x_4_x_T.y
CD3$x_4_x_x<-test2$x_4_x_x   + test2$x_4_x_x.x  + test2$x_4_x_x.y

CD3$x_x_13_T<-test2$x_x_13_T + test2$x_x_13_T.x + test2$x_x_13_T.y
CD3$x_x_13_x<-test2$x_x_13_x + test2$x_x_13_x.x + test2$x_x_13_x.y
CD3$x_x_x_T<-test2$x_x_x_T   + test2$x_x_x_T.x  + test2$x_x_x_T.y
CD3$x_x_x_x<-test2$x_x_x_x   + test2$x_x_x_x.x  + test2$x_x_x_x.y
CD3_compass<-CD3
CD3_compass$celltype<-"CD3";


dfs<-rbindlist(list(CD4_compass, CD8_compass, GD_compass), fill=TRUE)
dfs$TB<-factor(dfs$TB, levels=c("N", "HC", "LTBI", "TB"))
dfs<-arrange(dfs, celltype, TB)
###New attempts at doing it at once
dfs$unique<-paste(dfs$Donor, dfs$celltype)
metadata<-dplyr::select(dfs, Sample, Donor, TB, SM, Stim, celltype, unique)
metadata$group<-paste(metadata$TB,metadata$SM)

counts<- dplyr::select(dfs, G_4_13_T, 
                       G_4_x_T, G_x_13_T, G_4_13_x, x_4_13_T,
                       G_x_x_T,G_4_x_x,G_x_13_x, x_4_x_T, x_x_13_T, x_4_13_x,
                       G_x_x_x, x_x_x_T, x_4_x_x, x_x_13_x, 
                       x_x_x_x)
compass_nms<- c("IFNg&TNFa&IL4&IL13",
                    "IFNg&TNFa&IL4&!IL13", "IFNg&TNFa&!IL4&IL13", "IFNg&!TNFa&IL4&IL13", "!IFNg&TNFa&IL4&IL13",
                    "IFNg&TNFa&!IL4&!IL13","IFNg&!TNFa&IL4&!IL13","IFNg&!TNFa&!IL4&IL13","!IFNg&TNFa&IL4&!IL13","!IFNg&TNFa&!IL4&IL13", "!IFNg&!TNFa&IL4&IL13",
                    "IFNg&!TNFa&!IL4&!IL13", "!IFNg&TNFa&!IL4&!IL13", "!IFNg&!TNFa&IL4&!IL13", "!IFNg&!TNFa&!IL4&IL13",
                "!IFNg&!TNFa&!IL4&!IL13")

colnames(counts)<-compass_nms;

#SEA
inter<-intersect(metadata$unique[metadata$Stim =="UN"],metadata$unique[metadata$Stim =="SEA"])
n_u<-subset(counts,metadata$Stim =="UN" & metadata$unique %in% inter)
rownames(n_u)<-metadata$unique[metadata$Stim =="UN" & metadata$unique %in% inter]
n_s<-subset(counts,metadata$Stim =="SEA" & metadata$unique %in% inter)
rownames(n_s)<-metadata$unique[metadata$Stim =="SEA" & metadata$unique %in% inter]
metadata_sea<-data.frame(metadata[(metadata$Stim =="SEA" & metadata$unique %in% inter),])
fit_sea<-COMPASS::SimpleCOMPASS (n_s=n_s, n_u=n_u, meta=metadata_sea, individual_id = "unique", iterations = 10000, 
                                 replications = 8, verbose = TRUE);
save(fit_sea, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/all_sea2.RData")
write.csv(scores(fit_sea), 
          "/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_compass_sea.csv")

#SWAP
inter<-intersect(metadata$unique[metadata$Stim =="UN"],metadata$unique[metadata$Stim =="SWAP"])
n_u<-subset(counts,metadata$Stim =="UN" & metadata$unique %in% inter)
rownames(n_u)<-metadata$unique[metadata$Stim =="UN" & metadata$unique %in% inter]
n_s<-subset(counts,metadata$Stim =="SWAP" & metadata$unique %in% inter)
rownames(n_s)<-metadata$unique[metadata$Stim =="SWAP" & metadata$unique %in% inter]
metadata_SWAP<-data.frame(metadata[(metadata$Stim =="SWAP" & metadata$unique %in% inter),])
fit_swap<-COMPASS::SimpleCOMPASS (n_s=n_s, n_u=n_u, meta=metadata_SWAP, individual_id = "unique", iterations = 10000, 
                                 replications = 8, verbose = TRUE);
save(fit_swap, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/all_SWAP2.RData")
write.csv(scores(fit_swap), 
          "/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_compass_swap.csv")

#PMA
inter<-intersect(metadata$unique[metadata$Stim =="UN"],metadata$unique[metadata$Stim =="PMA"])
n_u<-subset(counts,metadata$Stim =="UN" & metadata$unique %in% inter)
rownames(n_u)<-metadata$unique[metadata$Stim =="UN" & metadata$unique %in% inter]
n_s<-subset(counts,metadata$Stim =="PMA" & metadata$unique %in% inter)
rownames(n_s)<-metadata$unique[metadata$Stim =="PMA" & metadata$unique %in% inter]
metadata_PMA<-data.frame(metadata[(metadata$Stim =="PMA" & metadata$unique %in% inter),])
fit_pma<-COMPASS::SimpleCOMPASS (n_s=n_s, n_u=n_u, meta=metadata_PMA, individual_id = "unique", iterations = 10000, 
                                  replications = 8, verbose = TRUE);
save(fit_pma, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/all_pma2.RData")

write.csv(scores(fit_pma), 
          "/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_compass_pma.csv")




####                OLD VERSION             #######################

####I know this works
compass_data<-CD4_compass
metadata<-dplyr::select(compass_data, Sample, Donor, TB, SM, Stim, celltype)
metadata$group<-paste(metadata$TB,metadata$SM)
counts<- dplyr::select(compass_data, G_4_13_T, G_4_13_x, G_4_x_T, G_4_x_x,
                       G_x_13_T, G_x_13_x, G_x_x_T, G_x_x_x, 
                       x_4_13_T, x_4_13_x, x_4_x_T, x_4_x_x, 
                       x_x_13_T, x_x_13_x, x_x_x_T, x_x_x_x)
compass_nms<- c(
    "IFNg&IL4&IL13&TNFa",
    "IFNg&IL4&IL13&!TNFa",
    "IFNg&IL4&!IL13&TNFa",
    "IFNg&IL4&!IL13&!TNFa",
    "IFNg&!IL4&IL13&TNFa",
    "IFNg&!IL4&IL13&!TNFa",
    "IFNg&!IL4&!IL13&TNFa",
    "IFNg&!IL4&!IL13&!TNFa",
    "!IFNg&IL4&IL13&TNFa",
    "!IFNg&IL4&IL13&!TNFa",
    "!IFNg&IL4&!IL13&TNFa",
    "!IFNg&IL4&!IL13&!TNFa",
    "!IFNg&!IL4&IL13&TNFa",
    "!IFNg&!IL4&IL13&!TNFa",
    "!IFNg&!IL4&!IL13&TNFa",
    "!IFNg&!IL4&!IL13&!TNFa")

colnames(counts)<-compass_nms
metadata$unique_id <- metadata$Donor

#UN
n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")

#WCL
n_wcl <- subset(counts, metadata$Stim =="WCL");
rownames(n_wcl) <- subset(metadata$unique_id,metadata$Stim=="WCL")
metadata_wcl<-data.frame(subset(metadata,metadata$Stim=="WCL"))

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id, metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_wcl))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_wcl<-COMPASS::SimpleCOMPASS (n_s=n_wcl, n_u=n_u, meta=metadata_wcl, individual_id = "unique_id", 
                iterations = 10000, replications = 8, verbose = TRUE)

#PEP
n_pep <- subset(counts, metadata$Stim =="PEP")
rownames(n_pep) <- subset(metadata$unique_id,metadata$Stim=="PEP")
metadata_pep<-data.frame(subset(metadata,metadata$Stim=="PEP"))

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_pep))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_pep<-COMPASS::SimpleCOMPASS (n_s=n_pep, n_u=n_u, meta=metadata_pep, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE)

#PMA
n_pma <- subset(counts, metadata$Stim =="PMA")
rownames(n_pma) <- subset(metadata$unique_id,metadata$Stim=="PMA")
metadata_pma<-data.frame(subset(metadata,metadata$Stim=="PMA"))

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_pma))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_pma<-COMPASS::SimpleCOMPASS (n_s=n_pma, n_u=n_u, meta=metadata_pma, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE);


#SEA
n_sea <- subset(counts, metadata$Stim =="SEA")
rownames(n_sea) <- subset(metadata$unique_id,metadata$Stim=="SEA")
metadata_sea<-data.frame(subset(metadata,metadata$Stim=="SEA"))

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_sea))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_sea<-COMPASS::SimpleCOMPASS (n_s=n_sea, n_u=n_u, meta=metadata_sea, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE);

#SWAP
n_swap<- subset(counts, metadata$Stim =="SWAP")
rownames(n_swap) <- subset(metadata$unique_id,metadata$Stim=="SWAP")
metadata_swap<-data.frame(subset(metadata,metadata$Stim=="SWAP"))

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_swap))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_swap<-COMPASS::SimpleCOMPASS (n_s=n_swap, n_u=n_u, meta=metadata_swap, individual_id = "unique_id", iterations = 10000, 
                                  replications = 8, verbose = TRUE);

save(fit_pma, "/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd4-pma.RData")
save(fit_pep, "/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd4-pep.RData")
save(fit_wcl, "/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd4-wcl.RData")
save(fit_sea, "/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd4-sea.RData")
save(fit_swap, "/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd4-swap.RData")

list<-c(fit_pma,fit_wcl,fit_pep,fit_sea,fit_swap)

#scores
scores_pma<-COMPASS::scores(fit_pma)
scores_pma<-dplyr::select(scores_pma, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_wcl<-COMPASS::scores(fit_wcl)
scores_wcl<-dplyr::select(scores_wcl, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_pep<-COMPASS::scores(fit_pep)
scores_pep<-dplyr::select(scores_pep, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_sea<-COMPASS::scores(fit_sea)
scores_sea<-dplyr::select(scores_sea, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_swap<-COMPASS::scores(fit_swap)
scores_swap<-dplyr::select(scores_swap, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)

dfs<-list(scores_pma, scores_pep, scores_wcl, scores_sea, scores_swap)
scores_cd4<-rbindlist(dfs)

######CD8
compass_data<-CD8_compass
metadata<-dplyr::select(compass_data, Sample, Donor, TB, SM, Stim, celltype)
metadata$group<-paste(metadata$TB,metadata$SM)
counts<- dplyr::select(compass_data, G_4_13_T, G_4_13_x, G_4_x_T, G_4_x_x,
                       G_x_13_T, G_x_13_x, G_x_x_T, G_x_x_x, 
                       x_4_13_T, x_4_13_x, x_4_x_T, x_4_x_x, 
                       x_x_13_T, x_x_13_x, x_x_x_T, x_x_x_x)
compass_nms<- c(
    "IFNg&IL4&IL13&TNFa",
    "IFNg&IL4&IL13&!TNFa",
    "IFNg&IL4&!IL13&TNFa",
    "IFNg&IL4&!IL13&!TNFa",
    "IFNg&!IL4&IL13&TNFa",
    "IFNg&!IL4&IL13&!TNFa",
    "IFNg&!IL4&!IL13&TNFa",
    "IFNg&!IL4&!IL13&!TNFa",
    "!IFNg&IL4&IL13&TNFa",
    "!IFNg&IL4&IL13&!TNFa",
    "!IFNg&IL4&!IL13&TNFa",
    "!IFNg&IL4&!IL13&!TNFa",
    "!IFNg&!IL4&IL13&TNFa",
    "!IFNg&!IL4&IL13&!TNFa",
    "!IFNg&!IL4&!IL13&TNFa",
    "!IFNg&!IL4&!IL13&!TNFa")

colnames(counts)<-compass_nms
metadata$unique_id <- metadata$Donor

#UN
n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")

#WCL
n_wcl <- subset(counts, metadata$Stim =="WCL");
rownames(n_wcl) <- subset(metadata$unique_id,metadata$Stim=="WCL")
metadata_wcl<-subset(metadata,metadata$Stim=="WCL")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id, metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_wcl))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_wcl<-COMPASS::SimpleCOMPASS (n_s=n_wcl, n_u=n_u, meta=metadata_wcl, individual_id = "unique_id", iterations = 10000, replications = 8, verbose = TRUE)

#PEP
n_pep <- subset(counts, metadata$Stim =="PEP")
rownames(n_pep) <- subset(metadata$unique_id,metadata$Stim=="PEP")
metadata_pep<-subset(metadata,metadata$Stim=="PEP")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_pep))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_pep<-COMPASS::SimpleCOMPASS (n_s=n_pep, n_u=n_u, meta=metadata_pep, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE)

#PMA
n_pma <- subset(counts, metadata$Stim =="PMA")
rownames(n_pma) <- subset(metadata$unique_id,metadata$Stim=="PMA")
metadata_pma<-subset(metadata,metadata$Stim=="PMA")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_pma))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_pma<-COMPASS::SimpleCOMPASS (n_s=n_pma, n_u=n_u, meta=metadata_pma, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE);


#SEA
n_sea <- subset(counts, metadata$Stim =="SEA")
rownames(n_sea) <- subset(metadata$unique_id,metadata$Stim=="SEA")
metadata_sea<-subset(metadata,metadata$Stim=="SEA")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_sea))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_sea<-COMPASS::SimpleCOMPASS (n_s=n_sea, n_u=n_u, meta=metadata_sea, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE);

#SWAP
n_swap<- subset(counts, metadata$Stim =="SWAP")
rownames(n_swap) <- subset(metadata$unique_id,metadata$Stim=="SWAP")
metadata_swap<-subset(metadata,metadata$Stim=="SWAP")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_swap))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_swap<-COMPASS::SimpleCOMPASS (n_s=n_swap, n_u=n_u, meta=metadata_swap, individual_id = "unique_id", iterations = 10000, 
                                  replications = 8, verbose = TRUE);

save(fit_pma, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd8-pma.RData")
save(fit_pep, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd8-pep.RData")
save(fit_wcl, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd8-wcl.RData")
save(fit_sea, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd8-sea.RData")
save(fit_swap, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/cd8-swap.RData")


list<-c(fit_pma,fit_wcl,fit_pep,fit_sea,fit_swap)

#scores
scores_pma<-COMPASS::scores(fit_pma)
scores_pma<-dplyr::select(scores_pma, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_wcl<-COMPASS::scores(fit_wcl)
scores_wcl<-dplyr::select(scores_wcl, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_pep<-COMPASS::scores(fit_pep)
scores_pep<-dplyr::select(scores_pep, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_sea<-COMPASS::scores(fit_sea)
scores_sea<-dplyr::select(scores_sea, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_swap<-COMPASS::scores(fit_swap)
scores_swap<-dplyr::select(scores_swap, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)

dfs<-list(scores_pma, scores_pep, scores_wcl, scores_sea, scores_swap)
scores_cd8<-rbindlist(dfs)



####GD 
compass_data<-GD_compass
metadata<-dplyr::select(compass_data, Sample, Donor, TB, SM, Stim, celltype)
metadata$group<-paste(metadata$TB,metadata$SM)
counts<- dplyr::select(compass_data, G_4_13_T, G_4_13_x, G_4_x_T, G_4_x_x,
                       G_x_13_T, G_x_13_x, G_x_x_T, G_x_x_x, 
                       x_4_13_T, x_4_13_x, x_4_x_T, x_4_x_x, 
                       x_x_13_T, x_x_13_x, x_x_x_T, x_x_x_x)
compass_nms<- c(
    "IFNg&IL4&IL13&TNFa",
    "IFNg&IL4&IL13&!TNFa",
    "IFNg&IL4&!IL13&TNFa",
    "IFNg&IL4&!IL13&!TNFa",
    "IFNg&!IL4&IL13&TNFa",
    "IFNg&!IL4&IL13&!TNFa",
    "IFNg&!IL4&!IL13&TNFa",
    "IFNg&!IL4&!IL13&!TNFa",
    "!IFNg&IL4&IL13&TNFa",
    "!IFNg&IL4&IL13&!TNFa",
    "!IFNg&IL4&!IL13&TNFa",
    "!IFNg&IL4&!IL13&!TNFa",
    "!IFNg&!IL4&IL13&TNFa",
    "!IFNg&!IL4&IL13&!TNFa",
    "!IFNg&!IL4&!IL13&TNFa",
    "!IFNg&!IL4&!IL13&!TNFa")

colnames(counts)<-compass_nms
metadata$unique_id <- metadata$Donor

#UN
n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")

#WCL
n_wcl <- subset(counts, metadata$Stim =="WCL");
rownames(n_wcl) <- subset(metadata$unique_id,metadata$Stim=="WCL")
metadata_wcl<-subset(metadata,metadata$Stim=="WCL")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id, metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_wcl))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_wcl<-COMPASS::SimpleCOMPASS (n_s=n_wcl, n_u=n_u, meta=metadata_wcl, individual_id = "unique_id", iterations = 10000, replications = 8, verbose = TRUE)

#PEP
n_pep <- subset(counts, metadata$Stim =="PEP")
rownames(n_pep) <- subset(metadata$unique_id,metadata$Stim=="PEP")
metadata_pep<-subset(metadata,metadata$Stim=="PEP")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_pep))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_pep<-COMPASS::SimpleCOMPASS (n_s=n_pep, n_u=n_u, meta=metadata_pep, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE)

#PMA
n_pma <- subset(counts, metadata$Stim =="PMA")
rownames(n_pma) <- subset(metadata$unique_id,metadata$Stim=="PMA")
metadata_pma<-subset(metadata,metadata$Stim=="PMA")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_pma))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_pma<-COMPASS::SimpleCOMPASS (n_s=n_pma, n_u=n_u, meta=metadata_pma, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE);


#SEA
n_sea <- subset(counts, metadata$Stim =="SEA")
rownames(n_sea) <- subset(metadata$unique_id,metadata$Stim=="SEA")
metadata_sea<-subset(metadata,metadata$Stim=="SEA")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_sea))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_sea<-COMPASS::SimpleCOMPASS (n_s=n_sea, n_u=n_u, meta=metadata_sea, individual_id = "unique_id", iterations = 10000, 
                                 replications = 8, verbose = TRUE);

#SWAP
n_swap<- subset(counts, metadata$Stim =="SWAP")
rownames(n_swap) <- subset(metadata$unique_id,metadata$Stim=="SWAP")
metadata_swap<-subset(metadata,metadata$Stim=="SWAP")

n_u<-subset(counts,metadata$Stim =="UN")
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN")
inter<-intersect(rownames(n_u),rownames(n_swap))
n_u<-subset(counts,metadata$Stim =="UN" & metadata$Donor %in% inter)
rownames(n_u) <- subset(metadata$unique_id,metadata$Stim=="UN" & metadata$Donor %in% inter)

fit_swap<-COMPASS::SimpleCOMPASS (n_s=n_swap, n_u=n_u, meta=metadata_swap, individual_id = "unique_id", iterations = 10000, 
                                  replications = 8, verbose = TRUE);

save(fit_pma, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/gd-pma.RData")
save(fit_pep, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/gd-pep.RData")
save(fit_wcl, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/gd-wcl.RData")
save(fit_sea, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/gd-sea.RData")
save(fit_swap, file="/Users/tarynam/Desktop/Flow-scripts/Computational Packages/COMPASS/gd-swap.RData")

list<-c(fit_pma,fit_wcl,fit_pep,fit_sea,fit_swap)

#scores
scores_pma<-COMPASS::scores(fit_pma)
scores_pma<-dplyr::select(scores_pma, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_wcl<-COMPASS::scores(fit_wcl)
scores_wcl<-dplyr::select(scores_wcl, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_pep<-COMPASS::scores(fit_pep)
scores_pep<-dplyr::select(scores_pep, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_sea<-COMPASS::scores(fit_sea)
scores_sea<-dplyr::select(scores_sea, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)
scores_swap<-COMPASS::scores(fit_swap)
scores_swap<-dplyr::select(scores_swap, Sample, Donor, TB, SM, Stim, celltype, FS, PFS)

dfs<-list(scores_pma, scores_pep, scores_wcl, scores_sea, scores_swap)
scores_gd<-rbindlist(dfs)

#GRAPHS
input<-list("WCL"=fit_wcl,"Peptide"=fit_pep)
plotCOMPASSResultStack(input,row_annotation = c("Stim","TB","SM") ,variable="Antigen")
plot(fit_wcl, order_by_max_functionality = FALSE, markers=c("IL13","IL4","IFNg","TNFa"), 
     row_annotation=c("TB","SM"), border_color=NA, fontsize=14, threshold=0)
#setting threshold to -1 will include everything
#the group variable isn't a factor but you can make it a factor to plot that group wise and do stats and stuff
#https://bioconductor.org/packages/3.7/bioc/manuals/COMPASS/man/COMPASS.pdf


dfs<-list(scores_cd4, scores_cd8, scores_gd)
compass_scores<-rbindlist(dfs, fill=TRUE)

write.csv(compass_scores, 
    "/Applications/Old Computer/Day Lab/Flow-Data/Clean/ICS_Compass_scores_complete_clean.csv")



