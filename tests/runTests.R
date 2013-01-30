

#' This is a modified version of BiocGenerics:::testPackage, now able to 
#' accept also a pattern for functions.
testPackage <- function (pkgname, subdir = "unitTests", 
                         pattern = "^test_.*\\.R$",
                         patternFun = "^test_.*$") 
{
    .failure_details <- function(result) {
        res <- result[[1L]]
        if (res$nFail > 0 || res$nErr > 0) {
            Filter(function(x) length(x) > 0, lapply(
                res$sourceFileResults, 
                function(fileRes) {
                    names(Filter(function(x) x$kind != "success", 
                                 fileRes))
                }))
        }
        else list()
    }
    require(pkgname, quietly = TRUE, character.only = TRUE) || 
        stop("package '", pkgname, "' not found")
    dir <- system.file(subdir, package = pkgname)
    if (nchar(dir) == 0L) 
        stop("unable to find unit tests, no '", subdir, "' dir")
    require("RUnit", quietly = TRUE) || stop("RUnit package not found")
    RUnit_opts <- getOption("RUnit", list())
    RUnit_opts$verbose <- 0L
    RUnit_opts$silent <- TRUE
    RUnit_opts$verbose_fail_msg <- TRUE
    options(RUnit = RUnit_opts)
    suite <- defineTestSuite(name = paste(pkgname, "RUnit Tests"), 
                             dirs = dir,
                             testFileRegexp = pattern, 
                             testFuncRegexp = patternFun,
                             rngKind = "default", 
                             rngNormalKind = "default")
    result <- runTestSuite(suite)
    cat("\n\n")
    printTextProtocol(result, showDetails = FALSE)
    if (length(details <- .failure_details(result)) > 0) {
        cat("\nTest files with failing tests\n")
        for (i in seq_along(details)) {
            cat("\n  ", basename(names(details)[[i]]), "\n")
            for (j in seq_along(details[[i]])) {
                cat("    ", details[[i]][[j]], "\n")
            }
        }
        cat("\n\n")
        stop("unit tests failed for package ", pkgname)
    }
    result
}

library(RUnit)

testPackage(pkgname="Prots4Prots", subdir = "unitTests", 
            pattern="^test_nrm.+\\.R",
            patternFun="^test_fun.+")
testPackage(pkgname="Prots4Prots", subdir = "unitTests", 
            pattern="^test_stt.+\\.R",
            patternFun="^test_fun.+")
testPackage(pkgname="Prots4Prots", subdir = "unitTests", 
            pattern="^test_mtc.+\\.R",
            patternFun="^test_fun.+")

sessionInfo()
