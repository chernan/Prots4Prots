## See Bioconductor guidelines for Unit testing
## http://www.bioconductor.org/developers/unitTesting-guidelines/

library(RUnit)

## -----------------------------------------------------------------------------
## Manually delete content of temp/ when running test suites of only one type
## -----------------------------------------------------------------------------

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

## Test graphic functions

testSuiteSTTGra <- 
    defineTestSuite("Statistic methods (STT functions)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_stt.+\\.R",
                    testFuncRegexp = "^test_gra.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")
testSuiteMTCGra <- 
    defineTestSuite("Multiple testing correction (MTC graphics)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_mtc.+\\.R",
                    testFuncRegexp = "^test_gra.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")

## Run all test suites

set.seed(0)

rUnitGraTestsOutput <- file.path(getwd(), 'temp', "runit_gra_output.pdf")
pdf(rUnitGraTestsOutput)
## Warning : 'testSuite' object is not of class 'RUnitTestSuite' is a bug
testResultGra <- runTestSuite(
    testSuites = list(
        STT = testSuiteSTTGra,
        MTC = testSuiteMTCGra
    ),
    verbose=0)
dev.off()

## View pdf output

file.show(rUnitGraTestsOutput)

## Report tests as HTML

rUnitGraTestsReportFile <- file.path(getwd(), 'temp', "runit_gra_tests_output.html")
printHTMLProtocol(testResultGra, fileName=rUnitGraTestsReportFile)
browseURL(paste0('file://', rUnitGraTestsReportFile) )

