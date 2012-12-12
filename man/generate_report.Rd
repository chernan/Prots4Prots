\name{generate_report}
\alias{generate_report}
\title{
Generate HTML report.
}
\description{
Generate a HTML report given a Rmd file. Intermediate md file is also created.
}
\usage{
generate_report(all.methods.label, output_file_name_Rmd, output_folder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all.methods.label}{
Label for the name of the output.
}
  \item{output_file_name_Rmd}{
Input Rmd file.
}
  \item{output_folder}{
Folder where output files will be created.
}
}
\details{
Multiple calls to knitr API.
}
\value{
Returns names of the HTML file.
}
\references{
http://yihui.name/knitr/
}
\author{
Celine Hernandez <Celine.Hernandez@unil.ch>
}
\section{Warning}{
Also creates duplicated intermediary md files.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (all.methods.label, output_file_name_Rmd, output_folder) 
{
    output_file_name_md <- paste(c(output_folder, "/output_", 
        all.methods.label, ".md"), collapse = "")
    output_file_name_html <- paste(c(output_folder, "/output_", 
        all.methods.label, ".html"), collapse = "")
    knit(output_file_name_Rmd, out = output_file_name_md)
    knitr::knit2html(output_file_name_md, output = output_file_name_html)
    return(output_file_name_html)
  }
}
\keyword{ knitr }
