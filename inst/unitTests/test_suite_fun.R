## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)

## -----------------------------------------------------------------------------
## Manually delete content of temp/ when running test suites of only one type
## -----------------------------------------------------------------------------

## Multiple Testing Correction

sapply(
    list.files(file.path(getwd(), 'R'), pattern="mtc_.+\\.R", full.names=TRUE),
    source
)

## Test numeric functions

testSuiteMTCFun <- 
    defineTestSuite("Multiple testing correction (MTC functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_mtc.+\\.R",
                    testFuncRegexp = "^test_fun.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")

## Run all test suites
## Warning : 'testSuite' object is not of class 'RUnitTestSuite' is a bug

testResultFun <- runTestSuite(
    testSuites = list(
        testSuiteMTCFun
    ), 
    verbose=0)

## Report tests as HTML

rUnitFunTestsReportFile <- file.path(getwd(), 'temp', "runit_fun_tests_output.html")
printHTMLProtocol(testResultFun, fileName=rUnitFunTestsReportFile)
browseURL(paste0('file://', rUnitFunTestsReportFile) )
