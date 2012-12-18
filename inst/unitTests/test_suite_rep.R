## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)
library(knitr)

## -----------------------------------------------------------------------------
## Manually delete content of temp/ when running test suites of only one type
## -----------------------------------------------------------------------------

## Multiple Testing Correction

sapply(
    list.files(file.path(getwd(), 'R'), pattern="mtc_.+\\.R", full.names=TRUE),
    source
)

## Test report functions

testSuiteMTCRep <- 
    defineTestSuite("Multiple testing correction (MTC report)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_mtc.+\\.R",
                    testFuncRegexp = "^test_rep.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")


## Run all test suites

rUnitRepTestsOutput <- file.path(getwd(), 'temp', "runit_rep_output.Rmd")
sink(rUnitRepTestsOutput, append=TRUE)
## Warning : 'testSuite' object is not of class 'RUnitTestSuite' is a bug
testResultRep <- runTestSuite(
    testSuites = list(
        MTC = testSuiteMTCRep
    ),
    verbose=0)
sink(file=NULL)

## Generate and view knitr report

rUnitRepTestsOutputMd <- 
    file.path(getwd(), 'temp', "runit_rep_output.md")
rUnitRepTestsOutputHtml <- 
    file.path(getwd(), 'temp', "runit_rep_output.html")

knit(rUnitRepTestsOutput, output=rUnitRepTestsOutputMd)
knitr::knit2html(rUnitRepTestsOutputMd, output=rUnitRepTestsOutputHtml)
browseURL(paste0('file://', rUnitRepTestsOutputHtml) )

## Report tests as HTML

rUnitRepTestsReportFile <- file.path(getwd(), 'temp', "runit_rep_tests_output.html")
printHTMLProtocol(testResultRep, fileName=rUnitRepTestsReportFile)
browseURL(paste0('file://', rUnitRepTestsReportFile) )

