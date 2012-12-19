
test_rep_generateKnitrReport <- function() {
    
    ## Setup dataset
    rmdFile <- file.path(getwd(), 'temp', 'output_test_rep_generateKnitrReport.Rmd')
    file.create(rmdFile)
    cat(
        'Test for R Markdown',
        '========================================================',
        '',
        'This is an R Markdown document, inspired from the default R Markown document created by Rstudio.', 
        'Markdown is a simple formatting syntax for authoring web pages.',
        '',
        'Kniting a R Markdown document generates a HTML page that includes both text content as well as the output of any embedded R code chunks within the document.', 
        'You can embed an R code chunk like this:',
        '',
        '```{r}',
        'summary(cars)',
        '```',
        '',
        'You can also embed plots, for example:',
        '',
        '```{r fig.width=7, fig.height=6}',
        'plot(cars)',
        '```',
        sep = "\n",
        file=rmdFile
    )
    
    
    ## Run function
    outputFileNameHtml <- generateKnitrReport(
        'test_rep_generateKnitrReport', 
        rmdFile, 
        file.path(getwd(), 'temp'))
    
    ## Test characteristics of output file
    checkTrue(file.exists(outputFileNameHtml))
    checkTrue(file.info(outputFileNameHtml)$size != 0)
}


