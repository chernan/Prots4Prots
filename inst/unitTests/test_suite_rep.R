## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)
library(knitr)
library(markdown)

## -----------------------------------------------------------------------------
## Manually delete content of temp/ when running test suites of only one type
## -----------------------------------------------------------------------------

## Normalization methods
sapply(
    list.files(file.path(getwd(), 'R'), pattern="nrm_.+\\.R", full.names=TRUE),
    source
)

## Statistic methods
sapply(
    list.files(file.path(getwd(), 'R'), pattern="stt_.+\\.R", full.names=TRUE),
    source
)

## Multiple Testing Correction
sapply(
    list.files(file.path(getwd(), 'R'), pattern="mtc_.+\\.R", full.names=TRUE),
    source
)

## Report generation with knitr
sapply(
    list.files(file.path(getwd(), 'R'), pattern="rep_.+\\.R", full.names=TRUE),
    source
)

## Test report functions

testSuiteNRMRep <- 
    defineTestSuite("Statistic methods (NRM functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_nrm.+\\.R",
                    testFuncRegexp = "^test_rep.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")
testSuiteSTTRep <- 
    defineTestSuite("Statistic methods (STT functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_stt.+\\.R",
                    testFuncRegexp = "^test_rep.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")
testSuiteMTCRep <- 
    defineTestSuite("Multiple testing correction (MTC report)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_mtc.+\\.R",
                    testFuncRegexp = "^test_rep.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")
testSuiteREPRep <- 
    defineTestSuite("Report generation with knitr (REP report)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_rep.+\\.R",
                    testFuncRegexp = "^test_rep.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")

## Reset knitr options to default
opts_knit$restore()

## Set needed option values
##  progress If knitr sould show progress bar. Set to FALSE because of sink
##   redirection during execution of test suites.
##  base.dir Location of figure folder. Needed.
opts_knit$set(progress=FALSE,
              base.dir=file.path(getwd(), 'temp')
)

## Run all test suites

set.seed(0)

rUnitRepTestsOutput <- file.path(getwd(), 'temp', "runit_rep_output.Rmd")
sink(rUnitRepTestsOutput)
## Warning : 'testSuite' object is not of class 'RUnitTestSuite'. Bug?
# suppressWarnings(
testResultRep <- runTestSuite(
    testSuites = list(
        NRM = testSuiteNRMRep,
        STT = testSuiteSTTRep,
        MTC = testSuiteMTCRep,
        REP = testSuiteREPRep
    ),
    verbose=0)
# )
sink(file=NULL)

## Generate and view knitr report

rUnitRepTestsOutputMd <- 
    file.path(getwd(), 'temp', "runit_rep_output.md")
knit(rUnitRepTestsOutput, output=rUnitRepTestsOutputMd)

rUnitRepTestsOutputHtml <- 
    file.path(getwd(), 'temp', "runit_rep_output.html")
markdownToHTML(file=rUnitRepTestsOutputMd, output=rUnitRepTestsOutputHtml)
browseURL(paste0('file://', rUnitRepTestsOutputHtml) )

## Reset knitr options to default
opts_knit$restore()

## Report tests as HTML

rUnitRepTestsReportFile <- file.path(getwd(), 'temp', "runit_rep_tests_output.html")
printHTMLProtocol(testResultRep, fileName=rUnitRepTestsReportFile)
browseURL(paste0('file://', rUnitRepTestsReportFile) )

