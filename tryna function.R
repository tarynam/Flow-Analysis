splitter<-function(data,column,splitter,one=NULL,two=NULL,three=NULL,four=NULL){
        library(data.table)
        
        labels<-tstrsplit(data[,column], splitter)
        data[,one]<-labels[[1]]
        data[,two]<-labels[[2]]
        data[,three]<-labels[[3]]
        data[,four]<-labels[[4]]
        return(data)
}

base.plot <- function(DF, group1, group2, value) {
        p <- ggplot(DF, aes(x=group1, y=value, col=group1))
        p <- p + theme_bw()
        p <- p + theme(legend.position="bottom")
        p <- p + geom_jitter(width=.1,height=0, shape=1,size=2)
        p <- p + scale_color_brewer(palette="Set1", guide=guide_legend(ncol=6, title=NULL))
        p <- p + xlab("") + ylab("")
        return(p)
}

replace.btwn.tables2<-function(datatable1, datatable2, id.var, ref.var, new.var, new.vals){
    library(dplyr)
    inter <- intersect(datatable1[,id.var], datatable2[,ref.var])
    tbl1a <- subset(datatable1, datatable1[,id.var] %in% inter)
    tbl1b <- subset(datatable1, !(datatable1[,id.var] %in% inter))
    tbl2  <- subset(datatable2, datatable2[,ref.var] %in% inter)
    tbl1a[,new.var] <- replace(x = tbl1a[,new.var],
                           list = (tbl1a[,id.var] %in% tbl2[,ref.var]),
                           values = as.character(tbl2[,new.vals]))
    new.table<-rbind(tbl1a, tbl1b)
}

pvals<-function(data, column){
    
    data<-filter(data, SM!="N")
    data$SM<-factor(data$SM, levels=c("X","SM"))
    data$TB<-factor(data$TB, levels=c("HC","LTBI","TB"))
    
    x<-split(data,data$TB) 
    A<-lapply(x, function(g) wilcox.test(g[,column]~g[,"SM"]))
    
    y<-split(data,data$SM) 
    B<-lapply(y, function(g) wilcox.test((g[,column]~g[,"TB"]), subset=g$TB %in% c("HC", "LTBI")))
    C<- lapply(y, function(g) wilcox.test((g[,column]~g[,"TB"]), subset=g$TB %in% c("LTBI", "TB")))
    D<- lapply(y, function(g) wilcox.test((g[,column]~g[,"TB"]), subset=g$TB %in% c("HC", "TB")))
    
    labels<-c("HC SM+ to HC SM-", "LTBI SM+ to LTBI SM-", "TB SM+ to TB SM-",
              "HC SM+ to LTBI SM+", "LTBI SM+ to TB SM+", "HC SM+ to TB SM+",
              "HC SM- to LTBI SM-", "LTBI SM- to TB SM-", "HC SM- to TB SM-"
    )
    pvals<-c(A$HC$p.value, A$LTBI$p.value, A$TB$p.value,
             B$SM$p.value, C$SM$p.value,  D$SM$p.value,
             B$X$p.value,  C$X$p.value,   D$X$p.value
    )
    table<-(cbind(labels,pvals))
    colnames(table)<-c("Comparison", paste(column,"p-values"," "))
    table
}

clean_flowjo_output<-function(datatable){
    library(dplyr)
    datatable$Sample<-gsub("1: |2: |3: |4: |5: |6: ", "", datatable$Sample)
    datatable$Sample<-gsub("XUN", "X_UN", datatable$Sample)
    datatable$Sample<-gsub("XPMA", "X_PMA", datatable$Sample)
    datatable$Sample<-gsub("XWCL", "X_WCL", datatable$Sample)
    datatable$Sample<-gsub("XPEP", "X_PEP", datatable$Sample)
    datatable$Sample<-gsub("XSEA", "X_SEA", datatable$Sample)
    datatable$Sample<-gsub("XSWAP", "X_SWAP", datatable$Sample)
    datatable$Sample<-gsub("_ ", "_", datatable$Sample)
    datatable$Sample<-gsub(" _", "_", datatable$Sample)
    datatable<-filter(datatable,!(Sample=="Mean"))%>%
        filter(!Sample=="StdDev")%>%
        filter(!Sample=="Sample")%>%
        filter(!Sample=="Sample")
    datatable[datatable=="¥"]<-NA
    datatable
    
}

splitter<-function(data,column,splitter,one=NULL,two=NULL,three=NULL,four=NULL){
    library(data.table)
    
    labels<-tstrsplit(data[,column], splitter)
    data[,one]<-labels[[1]]
    data[,two]<-labels[[2]]
    data[,three]<-labels[[3]]
    data[,four]<-labels[[4]]
    return(data)
}
#demo: HEL<-splitter(HEL, "StudyID", "-", one="StudyID", two="Visit", three="Iteration")

remove_samples<-function(datatable){
    library(dplyr)
    excluded_samples<-c("NK2126_HC_X_SEA_005.fcs", 
                        "NK2450_HC_SM_SEA_047.fcs" , 
                        "NK2402_LTBI_SM_SWAP_024.fcs",
                        "NK2342_TB_X_SWAP_006.fcs",
                        "NK2136_TB_X_WCL_003.fcs",
                        "NK2168_HC_SM_WCL_040.fcs",
                        "NK2421_LTBI_SM_PEP_010.fcs",
                        "NK2063_TB_X_SEA_005.fcs",
                        "NK2115_LTBI_X_PEP_010.fcs",
                        "NK2115_LTBI_X_WCL_009.fcs",
                        "NK2171_HC_X_WCL_009.fcs",
                        "NK2162_HC_X_SEA_023.fcs")
    datatable<-filter(datatable, !(Sample %in% excluded_samples))
    datatable
}

