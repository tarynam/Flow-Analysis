data(ICS)
E<-ConstructMIMOSAExpressionSet(ICS,
                                reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
                                other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
                                default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
                                .variables=.(TCELLSUBSET,CYTOKINE,UID),
                                featureCols=1,
                                ref.append.replace='_REF')

my.fitMCMC<-function(data, inits = NULL, iter = 250000, burn = 50000, thin = 1,
         tune = 100, outfile = basename(tempfile(tmpdir = ".", fileext = ".dat")),
         alternative = "greater", UPPER = 0.5, LOWER = 0.15, FAST = TRUE,
         EXPRATE = .1, pXi = c(1, 1), seed = 10){}

fixInNamespace(".fitMCMC", my.fitMCMC, ns="MIMOSA")

result2<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
               data=E, 
               method='EM',
               subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
               ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference', 
               RT=FALSE
)

