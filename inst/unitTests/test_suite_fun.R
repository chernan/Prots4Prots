## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)

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

## Test numeric functions

testSuiteNRMFun <- 
    defineTestSuite("Statistic methods (NRM functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_nrm.+\\.R",
                    testFuncRegexp = "^test_fun.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")
testSuiteSTTFun <- 
    defineTestSuite("Statistic methods (STT functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_stt.+\\.R",
                    testFuncRegexp = "^test_fun.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")
testSuiteMTCFun <- 
    defineTestSuite("Multiple testing correction (MTC functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_mtc.+\\.R",
                    testFuncRegexp = "^test_fun.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")

## Run all test suites
## Warning : 'testSuite' object is not of class 'RUnitTestSuite' is a bug

set.seed(0)

testResultFun <- runTestSuite(
    testSuites = list(
        NRM = testSuiteNRMFun,
        STT = testSuiteSTTFun,
        MTC = testSuiteMTCFun
    ), 
    verbose=0)

## Report tests as HTML

rUnitFunTestsReportFile <- file.path(getwd(), 'temp', "runit_fun_tests_output.html")
printHTMLProtocol(testResultFun, fileName=rUnitFunTestsReportFile)
browseURL(paste0('file://', rUnitFunTestsReportFile) )

