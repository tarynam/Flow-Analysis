data<-Cytokine[Cytokine$Stim=="PMA"]
data$mask = 0
data$mask[data$TH2.Freq > 1.5] = 1
data$mask = factor(data$mask, levels = c('1', '0'))

library(ggplot2)
base.plot <- function(DF, group1, group2, value) {
        p <- ggplot(DF, aes(x=group1, y=value, col=group1))
        p <- p + theme_bw()
        p <- p + theme(legend.position="bottom")
        p <- p + geom_jitter(width=.1,height=0, shape=1,size=2)
        p <- p + scale_color_brewer(palette="Set1", guide=guide_legend(ncol=6, title=NULL))
        p <- p + xlab("") + ylab("")
        return(p)
}

p <- base.plot(DF=data, SM, TB, Value)

p <- p + facet_grid(mask~TB, scales="free_y")
p <- p + theme(strip.text.y = element_blank())
#The faceting in the Y is backwards and I don't know why or how to change it
#I have tried releveling the mask variable and that didn't work
#And why is the bottom chunk huge compared to the top?
#I tried using the shrink and the space options in facet_grid 
#to make the bottom box smaller but it didn't work
