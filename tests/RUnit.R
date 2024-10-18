if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("gmresls", quietly=TRUE)) {
   library(RUnit)
   library(gmresls)

   testSuite <- defineTestSuite(
      name = "gmresls unit tests",
      dirs = system.file("unitTests", package = "gmresls"),
      testFuncRegexp = "^[Tt]est.+",
      rngKind = RNGkind()[1L]
   )
   Sys.setenv("R_TESTS"="")
   tests <- runTestSuite(testSuite)

   printTextProtocol(tests)

   if (getErrors(tests)$nFail > 0) stop("RUnit test failure")
   if (getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
