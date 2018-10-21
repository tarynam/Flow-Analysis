# another option
makePairs <- function(data) 
{
        grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
        grid <- subset(grid, x != y)
        all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
                xcol <- grid[i, "x"]
                ycol <- grid[i, "y"]
                data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
                           x = data[, xcol], y = data[, ycol], data)
        }))
        all$xvar <- factor(all$xvar, levels = names(data))
        all$yvar <- factor(all$yvar, levels = names(data))
        densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
                data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
        }))
        list(all=all, densities=densities)
}

# expand iris data frame for pairs plot
gg1 = makePairs(CD4[,3:18])

# new data frame mega iris
mega_CD4 = data.frame(gg1$all, Species=rep(CD4$TB, length=nrow(gg1$all)))

# pairs plot
ggplot(mega_CD4, aes_string(x = "x", y = "y")) + 
        facet_grid(xvar ~ yvar, scales = "free") + 
        geom_point(aes(colour=Species), na.rm = TRUE, alpha=0.8) + 
        stat_density(aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), 
                     data = gg1$densities, position = "identity", 
                     colour = "grey20", geom = "line")