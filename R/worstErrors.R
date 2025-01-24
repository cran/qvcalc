#' Accuracy of a Quasi-variance Approximation
#'
#' Computes the worst relative error, among all contrasts, for the standard
#' error as derived from a set of quasi variances.  For details of the method
#' see Menezes (1999) or Firth and Menezes (2004).
#'
#' @param qv.object An object of class \code{qv}
#'
#' @return A numeric vector of length 2, the worst negative relative error and
#' the worst positive relative error.
#'
#' @author David Firth, \email{d.firth@@warwick.ac.uk}
#'
#' @seealso \code{\link{qvcalc}}
#'
#' @references
#'
#' Firth, D. and Mezezes, R. X. de (2004) Quasi-variances.
#' \emph{Biometrika} \bold{91}, 69--80.
#' \doi{10.1093/biomet/91.1.65}
#'
#' McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
#' London: Chapman and Hall.
#'
#' Menezes, R. X. (1999) More useful standard errors for group and factor
#' effects in generalized linear models.  \emph{D.Phil. Thesis}, Department of
#' Statistics, University of Oxford.
#' 
#' @keywords regression models
#'
#' @examples
#' 
#' ##  Overdispersed Poisson loglinear model for ship damage data
#' ##  from McCullagh and Nelder (1989), Sec 6.3.2
#' library(MASS)
#' data(ships)
#' ships$year <- as.factor(ships$year)
#' ships$period <- as.factor(ships$period)
#' shipmodel <- glm(formula = incidents ~ type + year + period,
#'     family = quasipoisson,
#'     data = ships, subset = (service > 0), offset = log(service))
#' shiptype.qvs <- qvcalc(shipmodel, "type")
#' summary(shiptype.qvs, digits = 4)
#' worstErrors(shiptype.qvs)
#'
#' @export
#' 
worstErrors <- function(qv.object)
{
    reducedForm <- function(covmat, qvmat){
        nlevels <- dim(covmat)[1]
 	firstRow <- covmat[1, ]
 	ones <- rep(1, nlevels)
 	J <- outer(ones, ones)
 	notzero <- 2:nlevels
 	r.covmat <- covmat + (firstRow[1]*J) -
            outer(firstRow, ones) -
            outer(ones, firstRow)
 	r.covmat <- r.covmat[notzero, notzero]
 	qv1 <- qvmat[1, 1]
 	r.qvmat <- (qvmat + qv1*J)[notzero, notzero]
 	list(r.covmat, r.qvmat)}
    covmat <- qv.object$covmat
    qvmat <- diag(qv.object$qvframe$quasiVar)
    r.form <- reducedForm(covmat, qvmat)
    r.covmat <- r.form[[1]]
    r.qvmat <- r.form[[2]]
    inverse.sqrt <- solve(chol(r.covmat))
    evalues <- eigen(t(inverse.sqrt) %*% r.qvmat %*% inverse.sqrt,
                     symmetric=TRUE)$values
    sqrt(c(min(evalues), max(evalues))) - 1
}
