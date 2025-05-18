test_that("mzData works", {
  demoFile <- fs::path_package(
    "extdata",
    "WKL-0001.mzdata.xml",
    package = "readAgilent"
  )
  demoFile
  spectrData <- mzData$new(filename = demoFile)
  expect_equal(class(spectrData), c("mzData", "R6"))
  expect_equal(spectrData$length, 2)
  expect_equal(dim(spectrData$spectrum(number = 1)), c(25953, 2))
  expect_equal(dim(spectrData$spectrum(number = 2)), c(24947, 2))
  expect_equal(
    formatC(sum(spectrData$spectrum(number = 1)$mz, na.rm = TRUE), digits = 8),
    " 22070885"
  )
  expect_equal(
    formatC(sum(spectrData$spectrum(number = 2)$mz, na.rm = TRUE), digits = 8),
    " 22134568"
  )
  expect_equal(
    formatC(
      sum(spectrData$spectrum(number = 1)$intensity, na.rm = TRUE),
      digits = 8
    ),
    "  3006477"
  )
  expect_equal(
    formatC(
      sum(spectrData$spectrum(number = 2)$intensity, na.rm = TRUE),
      digits = 8
    ),
    "6564663.6"
  )
})
