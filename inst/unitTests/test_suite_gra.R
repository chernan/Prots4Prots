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

## Test graphic functions

testSuiteMTCGra <- 
    defineTestSuite("Multiple testing correction (MTC graphics)",
                    dirs = file.path(getwd(), 'inst', 'unitTests'),
                    testFileRegexp = "^test_mtc.+\\.R",
                    testFuncRegexp = "^test_gra.+",
                    rngKind = "Marsaglia-Multicarry",
                    rngNormalKind = "Kinderman-Ramage")

## Run all test suites

rUnitGraTestsOutput <- file.path(getwd(), 'temp', "runit_gra_output.pdf")
pdf(rUnitGraTestsOutput)
## Warning : 'testSuite' object is not of class 'RUnitTestSuite' is a bug
testResultGra <- runTestSuite(
    testSuites = list(
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

