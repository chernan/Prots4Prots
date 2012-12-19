library(knitr)

#' @title Generate HTML report with knitr.
#' 
#' @description
#' 
#' @details
#' Generate a HTML report given a Rmd file. Intermediate md file is also 
#' created.
#' TODO : use name of Rmd file instead of a label.
#' 
#' @references http://yihui.name/knitr/
#' @keywords knitr
#' 
#' @param outputLabel Label included in the name of the output file.
#' @param outputFileNameRmd Input Rmd file.
#' @param outputFolder Folder where output files will be created.
#' @return Full path of the generated HTML file.
generateKnitrReport <- function(outputLabel, outputFileNameRmd, outputFolder) {
    
    outputFileNameMd <- 
        paste(c(outputFolder, '/output_', outputLabel, '.md'), 
              collapse='')
    outputFileNameHtml <- 
        paste(c(outputFolder, '/output_', outputLabel, '.html'), 
              collapse='')
    
    ## Set needed option values
    ## See http://yihui.name/knitr/options#package_options
    ##  base.dir Location of figure folder. Needed.
    ##  out.format Automatically detected
    ##  root.dir Automatically detected
    ##  output.dir Automatically detected
    oldBasedir <- opts_knit$get("base.dir")
    opts_knit$set(base.dir=outputFolder)
    
    ## First step, transform R markdown to markdown
    knit(outputFileNameRmd, output=outputFileNameMd)
    ## Then only transform markdown to HTML
    markdownToHTML(file=outputFileNameMd, output=outputFileNameHtml)

    ## Reset option values
    opts_knit$set(base.dir=oldBasedir)
    
    return(outputFileNameHtml)
}