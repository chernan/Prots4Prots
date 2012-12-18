## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)

## Setup and clean temp folder for outputs

dir.create(file.path(getwd(), '/temp'), showWarnings=FALSE)
sapply(
    list.files(file.path(getwd(), '/temp'), pattern="\\.Rmd|\\.md|\\.txt|\\.html", full.names=TRUE),
    file.remove
)

## Multiple Testing Correction

sapply(
    list.files(file.path(getwd(), 'R'), pattern="mtc_.+\\.R", full.names=TRUE),
    source
)
testsuite.mtc <- 
    defineTestSuite("mtc",
                     dirs = file.path(getwd(), 'inst', 'unitTests'),
                     testFileRegexp = "^test_mtc.+\\.R",
                     testFuncRegexp = "^test_.+",
                     rngKind = "Marsaglia-Multicarry",
                     rngNormalKind = "Kinderman-Ramage")
testResult <- runTestSuite(testsuite.mtc, verbose=0)


## Report tests as HTML

testReportFile <- file.path(getwd(), 'temp', "test_output.html")
printHTMLProtocol(testResult, fileName=testReportFile)
browseURL(paste0('file://', testReportFile) )

