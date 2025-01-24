testthat::test_that("qvcalc.default works for glm", {
  ##  Overdispersed Poisson loglinear model for ship damage data
  ##  from McCullagh and Nelder (1989), Sec 6.3.2
  testthat::skip_if_not_installed("MASS")
  require(MASS)

  data(ships)
  ships$year <- as.factor(ships$year)
  ships$period <- as.factor(ships$period)
  shipmodel <- glm(formula = incidents ~ type + year + period,
                   family = quasipoisson,
                   data = ships,
                   subset = (service > 0),
                   offset = log(service))
  shiptype.qv <- qvcalc(shipmodel, "type")

  testthat::expect_named(
    shiptype.qv,
    c("covmat", "qvframe", "dispersion", "relerrs", "factorname",
      "coef.indices", "modelcall")
  )

  testthat::expect_equal(
    shiptype.qv$qvframe$quasiSE,
    c(0.200965853501369, 0.112705161723464,
      0.375250492322657, 0.323893536427565,
      0.232171862907782)
  )

  testthat::expect_error({
    qvcalc(shipmodel, "period")
  }, regexp = "factors with 3 or more levels")

  testthat::expect_error({
    qvcalc(shipmodel)
  }, regexp = "both NULL")

  summary(shiptype.qv)
  print(shiptype.qv)
  testthat::expect_silent({
    plot(shiptype.qv)
  })

  # no intercept model
  shipmodel2 <- glm(formula = incidents ~ type + year + period - 1,
                   family = quasipoisson,
                   data = ships,
                   subset = (service > 0),
                   offset = log(service))
  shiptype.qv2 <- qvcalc(shipmodel2, "type")

  testthat::expect_equal(
    shiptype.qv$qvframe$quasiSE,
    shiptype.qv2$qvframe$quasiSE
  )

  # no intercept model - using coef indices
  shiptype.qv2 <- qvcalc(shipmodel2, coef.indices = 1:5)

  testthat::expect_equal(
    shiptype.qv$qvframe$quasiSE,
    shiptype.qv2$qvframe$quasiSE
  )

})
