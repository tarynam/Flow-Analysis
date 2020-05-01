counts<-read.csv("/Applications/Old Computer/Day Lab/Flow-Data/clean/OG_counts_summary_clean.csv")

pkgs <- c("MIMOSA","Biobase","ggplot2","MASS", "plyr","reshape")
    # see what packages are currently installed
installed_pacakges <- row.names(installed.packages()) # loop over the needed packages
for(p in pkgs){
        # check if package is installed
        already_installed <- p %in% installed_pacakges
        # if not already installed, install it
        if(!already_installed){ install.packages(p)
        }
        # and load package
        library(p, character.only = TRUE)}    
    

#When I first tried to sub in a new value for EXPRATE, I couldn't find the actual function
#That .fitMCMC executes so I left the brackets empty
my.fitMCMC<-function (data, inits = NULL, iter = 250000, burn = 50000, thin = 1, 
          tune = 100, outfile = basename(tempfile(tmpdir = ".", fileext = ".dat")), 
          alternative = "greater", UPPER = 0.5, LOWER = 0.15, FAST = TRUE, 
          EXPRATE = .01, pXi = c(1, 1), seed = 10)
{
    set.seed(seed)
    alternative <- match.arg(alternative, c("greater", "not equal"))
    data <- MIMOSA:::icsdata2mvicsdata(data)
    if (is.null(inits)) {
        r <- MDMix(data)
        inits <- list(alpha.s = r@par.stim, alpha.u = r@par.unstim, 
                      q = r@w[1], z = round(r@z))
    }
    if (alternative == "greater") {
        ps <- t(do.call(cbind, apply(data$n.stim, 1, function(x) (data.frame(prop.table(x))[-1L, 
                                                                                            , drop = FALSE]))))
        pu <- t(do.call(cbind, apply(data$n.unstim, 1, function(x) (data.frame(prop.table(x))[-1L, 
                                                                                              , drop = FALSE]))))
        filter <- sapply(1:nrow(ps), function(i) all(ps[i, ] <= 
                                                         pu[i, ]))
        FILTER = TRUE
    }
    else {
        filter <- rep(FALSE, nrow(data$n.stim))
        FILTER <- FALSE
    }
    result <- .Call("C_fitMCMC", as.matrix(data$n.stim), as.matrix(data$n.unstim), 
                    as.vector(inits$alpha.s), as.vector(inits$alpha.u), as.vector(inits$q), 
                    as.matrix(inits$z), as.vector(iter), as.vector(burn), 
                    as.vector(thin), as.numeric(tune), as.character(outfile), 
                    as.vector(filter), as.numeric(UPPER), as.numeric(LOWER), 
                    FILTER, FAST, as.numeric(EXPRATE), as.numeric(pXi))
    if (inherits(result, "character")) {
        return(result)
    }
    result$z <- cbind(result$z, 1 - result$z)
    result$getmcmc <- function(x = outfile) {
        coda:::mcmc(read.table(x, sep = "\t", header = TRUE))
    }
    result$getP <- function(x = paste(outfile, "P", sep = ""), 
                            thin = 1) {
        if (thin > 1) {
            nc <- length(strsplit(readLines(x, 1), "\t")[[1]])
            thins <- paste("p", paste(rep(";n", thin - 1), collapse = ""), 
                           sep = "")
            s <- sprintf("sed -n '%s' %s|cut -f %s-%s", thins, 
                         x, (nc/3 + 1), nc)
            con <- pipe(s)
            d <- do.call(rbind, lapply(strsplit(readLines(con), 
                                                "\t")[-1L], as.numeric))
            colnames(d) <- strsplit(readLines(x, 1), "\t")[[1]][(nc/3 + 
                                                                     1):nc]
            d <- split(as.list(data.frame(d)), gl(nc/3, 2))
            close(con)
        }
        else {
            d <- coda:::mcmc(read.table(x, sep = "\t", header = TRUE))
            nc <- ncol(d)
            d <- split(as.list(data.frame(d[, (nc/3 + 1):nc])), 
                       gl(nc/3, 2))
        }
        d
    }
    attr(result, "class") <- c(attr(result, "class"), "MDMixResult")
    attr(result, "pData") <- attr(data, "pData")
    result$n.stim <- data$n.stim
    result$n.unstim <- data$n.unstim
    result
}

#The function reads in fine and I can both assign and fix the new function into MIMOSA's name space
assignInNamespace(".fitMCMC", my.fitMCMC, ns="MIMOSA")

## CD4
data<-dplyr::select(counts, Donor, TB, SM, Stim, CD4_count, CD4_neg)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim %in% 'UN',
                                measure.columns=c('CD4_neg','CD4_count'),
                                other.annotations=c('Stim', 'Donor', 'TB', 'SM'),
                                default.cast.formula=component~Donor+Stim,
                                .variables=.(Donor))

test <- MIMOSA(CD4_neg+CD4_count~Donor|Stim,
               data=E,
               method="mcmc")

CD4.probabilities<-MIMOSA::getZ(test)
CD4.counts<-MIMOSA::countsTable(test)
CD4.fdr<-MIMOSA::fdr(test)
CD4.Results<-data.frame(cbind(CD4.counts,CD4.probabilities,CD4.fdr))
CD4.Results$id<-row.names(CD4.counts)
CD4.Results$id<-gsub("_Treatment","", CD4.Results$id)
CD4.Results$celltype<-"CD4"
names(CD4.Results)<-gsub("CD4_", "", names(CD4.Results))
CD4.summary<-MIMOSA::getW(test)

## CD8
data<-dplyr::select(counts, Donor, TB, SM, Stim, CD8_count, CD8_neg)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim %in% 'UN',
                                measure.columns=c('CD8_neg','CD8_count'),
                                other.annotations=c('Stim', 'Donor', 'TB', 'SM'),
                                default.cast.formula=component~Donor+Stim,
                                .variables=.(Donor))

test <- MIMOSA(CD8_neg+CD8_count~Donor|Stim,
               data=E,
               method="mcmc")

CD8.probabilities<-MIMOSA::getZ(test)
CD8.counts<-MIMOSA::countsTable(test)
CD8.fdr<-MIMOSA::fdr(test)
CD8.Results<-data.frame(cbind(CD8.counts, CD8.probabilities, CD8.fdr))
CD8.Results$id<-row.names(CD8.counts)
CD8.Results$id<-gsub("_Treatment","", CD8.Results$id)
CD8.Results$celltype<-"CD8"
names(CD8.Results)<-gsub("CD8_", "", names(CD8.Results))
CD8.summary<-MIMOSA::getW(test)


## GD
data<-dplyr::select(counts, Donor, TB, SM, Stim, GD_count, GD_neg)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim %in% 'UN',
                                measure.columns=c('GD_neg','GD_count'),
                                other.annotations=c('Stim', 'Donor', 'TB', 'SM'),
                                default.cast.formula=component~Donor+Stim,
                                .variables=.(Donor))

test <- MIMOSA(GD_neg+GD_count~Donor|Stim,
               data=E,
               method="mcmc")

GD.probabilities<-MIMOSA::getZ(test)
GD.counts<-MIMOSA::countsTable(test)
GD.fdr<-MIMOSA::fdr(test)
GD.Results<-data.frame(cbind(GD.counts, GD.probabilities, GD.fdr))
GD.Results$id<-row.names(GD.counts)
GD.Results$id<-gsub("_Treatment","", GD.Results$id)
GD.Results$celltype<-"GD"
names(GD.Results)<-gsub("GD_", "", names(GD.Results))
GD.summary<-MIMOSA::getW(test)


## CD3
data<-dplyr::select(counts, Donor, TB, SM, Stim, CD3_count, CD3_neg)
E<-ConstructMIMOSAExpressionSet(data,
                                reference=Stim %in% 'UN',
                                measure.columns=c('CD3_neg','CD3_count'),
                                other.annotations=c('Stim', 'Donor', 'TB', 'SM'),
                                default.cast.formula=component~Donor+Stim,
                                .variables=.(Donor))

test <- MIMOSA(CD3_neg+CD3_count~Donor|Stim,
               data=E,
               method="mcmc")

CD3.probabilities<-MIMOSA::getZ(test)
CD3.counts<-MIMOSA::countsTable(test)
CD3.fdr<-MIMOSA::fdr(test)
CD3.Results<-data.frame(cbind(CD3.counts, CD3.probabilities, CD3.fdr))
CD3.Results$id<-row.names(CD3.counts)
CD3.Results$id<-gsub("_Treatment","", CD3.Results$id)
CD3.Results$celltype<-"CD3"
names(CD3.Results)<-gsub("CD3_", "", names(CD3.Results))
CD3.summary<-MIMOSA::getW(test)


MIMOSA.OG<-rbind(CD3.Results, CD4.Results, CD8.Results, GD.Results)

MIMOSA.OG$include<-0
for(i in 1:nrow(MIMOSA.OG)){
    if(MIMOSA.OG$Pr.response[i]>0.7 & MIMOSA.OG$fdr[i]<0.03){
        MIMOSA.OG$include[i]<-1
        print(MIMOSA.OG$include[i])}
}

MIMOSA.OG$include2<-0
for(i in 1:nrow(MIMOSA.OG)){
    if(MIMOSA.OG$Pr.response[i]>0.7 & MIMOSA.OG$fdr[i]<0.05){
        MIMOSA.OG$include2[i]<-1
        print(MIMOSA.OG$include2[i])}
}

write.csv(MIMOSA.OG, "/Applications/Old Computer/Day Lab/Flow-Data/clean/OG_MIMOSA_complete.csv")

