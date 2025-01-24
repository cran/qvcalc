testthat::test_that("qvcalc.itempar works for raschmodel", {
  testthat::skip_if_not_installed("psychotools")
  require(psychotools)

  data("VerbalAggression", package = "psychotools")
  raschmod <- raschmodel(VerbalAggression$resp2)
  ip1 <- itempar(raschmod)
  qv1 <- qvcalc(ip1)

  testthat::expect_named(
    qv1,
    c("covmat", "qvframe", "dispersion", "relerrs", "factorname",
      "coef.indices", "modelcall")
  )

  testthat::expect_equal(
    qv1$qvframe$quasiSE,
    c(0.142789774392438, 0.142789905559146, 0.132814807950912, 0.131468016260111,
      0.130384669951173, 0.1376080489105, 0.157119357966682, 0.136496922430792,
      0.134312458828237, 0.13043098568674, 0.130367526373492, 0.15154135165623,
      0.132499539340459, 0.130872793612835, 0.134902557606665, 0.152212972001133,
      0.152900324421324, 0.229780852197773, 0.139997658112739, 0.134312458828237,
      0.131624427111856, 0.131867589722712, 0.140748230935064, 0.170242788920416
    )
  )

  ip2 <- itempar(raschmod, alias = FALSE)
  testthat::expect_error({
      qvcalc(ip2)
  })


  ip2 <- itempar(raschmod, vcov = FALSE)
  testthat::expect_error({
    qvcalc(ip2)
  }, regexp = "vcov")


})
