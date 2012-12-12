library("ggplot2")
# library("MASS")

display_linearreg <- function(data.values, indexX=1, indexY=2, annX=10, annY=19, show00=FALSE) {
    model <- lm(data.values[,indexY]~data.values[,indexX])
    
    eq <- paste(
             c(
                 'y = ', 
                 format(coef(model)[1], digits=2),
                 ifelse( (coef(model)[2]>=0) , ' + ', ' - '),
                 format(abs(coef(model)[2]), digits=3),
                 ' . x'
             ), collapse=''
    )    
    r2val <- paste( c('r2 = ',
                      format(summary(model)$r.squared, digits=3)
    ), collapse='')    
    p <- ggplot(data=data.values, aes_string(x=names(data.values)[indexX], y=names(data.values)[indexY]))
    if(isTRUE(show00)) {
        p <- p + geom_point() + geom_hline(yintercept=0, colour="grey") + geom_vline(xintercept=0, colour="grey")
    }
    p + geom_point() + 
        geom_smooth(method=lm, formula=y~x, se=FALSE, color="blue") +
        annotate("text", x=annX, y=annY+0.2, hjust=0, vjust=0, label=eq) + 
        annotate("text", x=annX, y=annY-0.2, hjust=0, vjust=1, label=r2val) 
    
}

