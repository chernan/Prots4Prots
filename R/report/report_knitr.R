library(knitr)
library(markdown)

generate_report <- function(all.methods.label, output_file_name_Rmd, output_folder) {
    
    
    output_file_name_md <- paste(c(output_folder,'/output_',all.methods.label,'.md'), collapse='')
    output_file_name_html <- paste(c(output_folder,'/output_',all.methods.label,'.html'), collapse='')    
    
    knit(output_file_name_Rmd, out=output_file_name_md)
    knitr::knit2html(output_file_name_md, output=output_file_name_html)
    #     markdownToHTML(output_file_name_md, output=output_file_name_html)
    #     browseURL(output_file_name_html)
    
    return(output_file_name_html)
    
}