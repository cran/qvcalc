qvcalc.itempar <- function(object, ...){
  if (!(attr(object, "alias"))) stop(
                  "the itempar object was not built with 'alias = TRUE'")
  vc <- vcov(object)
  if (any(is.na(vc))) stop("the itempar object was not built with 'vcov = TRUE'")
  qvcalc.default(vc, estimates = object)
}
