
testthat::test_that("qvcalc.survival works", {
  testthat::skip_if_not_installed("survival")
  require(survival)
  ## Example of a "coxph" model
  cancer_model <- coxph(Surv(time,status) ~ celltype,
                        data = survival::veteran)
  celltype_qv <- qvcalc(cancer_model, "celltype")

  testthat::expect_named(
    celltype_qv,
    c("covmat", "qvframe", "dispersion", "relerrs", "factorname",
      "coef.indices", "modelcall")
  )

  testthat::expect_equal(
    celltype_qv$qvframe$quasiSE,
    c(0.202377154559161, 0.150182798946618,
      0.204942728570679, 0.199186083176635)
  )

})


## Example of a "survreg" model
testthat::test_that("qvcalc.survreg works", {
  testthat::skip_if_not_installed("survival")
  require(survival)
  cancer_model <- survreg(Surv(time,status) ~ celltype,
                           data = survival::veteran,
                           dist = "weibull")
  celltype_qv <- qvcalc(cancer_model, "celltype")
  testthat::expect_named(
    celltype_qv,
    c("covmat", "qvframe", "dispersion", "relerrs", "factorname",
      "coef.indices", "modelcall")
  )

  testthat::expect_equal(
    celltype_qv$qvframe$quasiSE,
    c(0.185265780005981, 0.153700051379706,
      0.202218015403209, 0.202439832203988)
  )
})
