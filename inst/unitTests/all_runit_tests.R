## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)

## -----------------------------------------------------------------------------
## NB:Manually delete content of temp/ when running test suites of only one type
## -----------------------------------------------------------------------------

## Setup and clean temp/ folder for outputs

dir.create(file.path(getwd(), '/temp'), showWarnings=FALSE)
sapply(
    list.files(
        file.path(getwd(), '/temp'), 
        pattern = "\\.Rmd|\\.md|\\.txt|\\.html|\\.pdf", 
        full.names = TRUE),
    file.remove
)

## Run all types of test suites: Fun(ctions), Gra(phics) and Rep(ort)

sapply(
    list.files(file.path(getwd(), 'inst', 'unitTests'), 
               pattern="test_suite_.+\\.R", full.names=TRUE),
    source
)

