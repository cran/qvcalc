##  The qvcalc generic and default method.  And methods for print, summary and plot,
##  for a resulting "qv" object.


#' Quasi Variances for Model Coefficients
#'
#' Computes a set of quasi variances (and corresponding quasi standard errors)
#' for estimated model coefficients relating to the levels of a categorical
#' (i.e., factor) explanatory variable.  For details of the method see Firth
#' (2000), Firth (2003) or Firth and de Menezes (2004).  Quasi variances
#' generalize and improve the accuracy of \dQuote{floating absolute risk}
#' (Easton et al., 1991).  This device for economical model summary was first
#' suggested by Ridout (1989).
#'
#'
#' The \code{qvcalc.default} method is the computational backend for all other,
#' class-specific methods.
#'
#' In \code{qvcalc.default}, none of the arguments other than \code{object} is
#' used in computing the result.  The remaining arguments are simply passed
#' through to the result object as components to help with record-keeping etc.
#'
#' In \code{qvcalc.lm}, at least one of \code{factorname} or
#' \code{coef.indices} must be non-\code{NULL}.  The value of
#' \code{coef.indices}, if non-\code{NULL}, determines which rows and columns
#' of the model's variance-covariance matrix to use.  If \code{coef.indices}
#' contains a zero, then an extra row and column are included at the indicated
#' position, to represent the zero variances and covariances associated with a
#' reference level.  If \code{coef.indices} is \code{NULL}, then
#' \code{factorname} should be the name of a factor effect in the model, and is
#' used in order to extract the necessary variance-covariance estimates.
#'
#' For \code{qvcalc.itempar}, the \code{"itempar"} object must have the full
#' variance-covariance matrix in its \code{"vcov"} attribute, and must have its
#' \code{"alias"} attribute be \code{TRUE}.  These attributes result from use
#' of the default arguments \code{vcov = TRUE, alias = TRUE} when the
#' \code{\link[psychotools]{itempar}} function is called.
#'
#' Ordinarily the quasi variances are positive and so their square roots (the
#' quasi standard errors) exist and can be used in plots, etc.
#'
#' Occasionally one (and only one) of the quasi variances is negative, and so
#' the corresponding quasi standard error does not exist (it appears as
#' \code{NaN}).  This is fairly rare in applications, and when it occurs it is
#' because the factor of interest is strongly correlated with one or more other
#' predictors in the model.  It is not an indication that quasi variances are
#' inaccurate.  An example is shown below using data from the \code{car}
#' package: the quasi variance approximation is exact (since \code{type} has
#' only 3 levels), and there is a negative quasi variance.  The quasi variances
#' remain perfectly valid (they can be used to obtain inference on any
#' contrast), but it makes no sense to plot `comparison intervals' in the usual
#' way since one of the quasi standard errors is not a real number.
#'
#' @aliases qvcalc qvcalc.default qvcalc.lm qvcalc.itempar qvcalc.coxph
#' qvcalc.survreg summary.qv print.qv
#' @param object For \code{qvcalc.default}, this is the covariance (sub)matrix
#' for the parameters of interest (including any that may have been constrained
#' to zero).  For the generic \code{qvcalc}, the \code{object} can be any
#' object for which the relevant S3 method has been defined.  These currently
#' include many types of regression model (via \code{qvcalc.lm}), including
#' objects of classes \code{\link[survival]{coxph}} and
#' \code{\link[survival]{survreg}}; and also objects of class
#' \code{\link[psychotools]{itempar}}.
#' @param factorname Either \code{NULL}, or a character vector of length 1
#' @param coef.indices Either \code{NULL}, or a numeric vector of length at
#' least 3
#' @param labels An optional vector of row names for the \code{qvframe}
#' component of the result (redundant if \code{object} is a model)
#' @param dispersion an optional scalar multiplier for the covariance matrix,
#' to cope with overdispersion for example
#' @param estimates an optional vector of estimated coefficients (redundant if
#' \code{object} is a model, for example)
#' @param modelcall optional, the call expression for the model of interest
#' (redundant if \code{object} is a model with its own \code{call} component)
#' @param ... other arguments to pass to \code{qv.default}
#' @return A list of class \code{qv}, with components \item{covmat}{the full
#' variance-covariance matrix for the estimated coefficients corresponding to
#' the factor of interest} \item{qvframe}{a data frame with variables
#' \code{estimate}, \code{SE}, \code{quasiSE} and \code{quasiVar}, the last two
#' being a quasi standard error and quasi-variance for each level of the factor
#' of interest} \item{relerrs}{relative errors for approximating the standard
#' errors of all simple contrasts} \item{factorname}{the factor name if given}
#' \item{coef.indices}{the coefficient indices if given} \item{modelcall}{if
#' \code{object} is a model, \code{object$call}; otherwise \code{NULL}}
#' @author David Firth, \email{d.firth@@warwick.ac.uk}
#' @seealso \code{\link{worstErrors}}, \code{\link{plot.qv}}.
#' @references Easton, D. F, Peto, J. and Babiker, A. G. A. G. (1991) Floating
#' absolute risk: an alternative to relative risk in survival and case-control
#' analysis avoiding an arbitrary reference group.  \emph{Statistics in
#' Medicine} \bold{10}, 1025--1035. \doi{10.1002/sim.4780100703}
#'
#' Firth, D. (2000) Quasi-variances in Xlisp-Stat and on the web.
#' \emph{Journal of Statistical Software} \bold{5.4}, 1--13.
#' \doi{10.18637/jss.v005.i04}
#'
#' Firth, D. (2003) Overcoming the reference category problem in the
#' presentation of statistical models. \emph{Sociological Methodology}
#' \bold{33}, 1--18. \doi{10.1111/j.0081-1750.2003.t01-1-00125.x}
#'
#' Firth, D. and de Mezezes, R. X. (2004) Quasi-variances.  \emph{Biometrika}
#' \bold{91}, 65--80. \doi{10.1093/biomet/91.1.65}
#'
#' McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
#' London: Chapman and Hall.
#'
#' Menezes, R. X. de (1999) More useful standard errors for group and factor
#' effects in generalized linear models.  \emph{D.Phil. Thesis}, Department of
#' Statistics, University of Oxford.
#'
#' Ridout, M.S. (1989). Summarizing the results of fitting generalized linear
#' models to data from designed experiments. In: \emph{Statistical Modelling:
#' Proceedings of GLIM89 and the 4th International Workshop on Statistical
#' Modelling held in Trento, Italy, July 17--21, 1989} (A. Decarli et al.,
#' eds.), pp 262--269. New York: Springer.
#' @keywords models regression
#' @examples
#'
#' ##  Overdispersed Poisson loglinear model for ship damage data
#' ##  from McCullagh and Nelder (1989), Sec 6.3.2
#' if (require(MASS)) {
#'     data(ships)
#'     ships$year <- as.factor(ships$year)
#'     ships$period <- as.factor(ships$period)
#'     shipmodel <- glm(formula = incidents ~ type + year + period,
#'                      family = quasipoisson,
#'                      data = ships,
#'                      subset = (service > 0),
#'                      offset = log(service))
#'     shiptype.qv <- qvcalc(shipmodel, "type")
#'
#'     ## We can plot "comparison intervals" as follows:
#'     ##   plot(shiptype.qv, xlab = "ship type")
#'
#'     ## An equivalent result by using the coef.indices argument instead:
#'     ##   shiptype.qv2 <- qvcalc(shipmodel, coef.indices = c(0, 2:5))
#'
#'     summary(shiptype.qv, digits = 4)
#' }
#'
#' ## Example of a "coxph" model
#' if(require(survival)) {
#'     data("veteran", package = "survival")
#'     cancer_model <- coxph(Surv(time,status) ~ celltype, data = veteran)
#'     celltype_qv <- qvcalc(cancer_model, "celltype")
#'     summary(celltype_qv)
#' }
#'
#' ## Example of a "survreg" model
#' if(require(survival)) {
#'     data("veteran", package = "survival")
#'     cancer_model2 <- survreg(Surv(time,status) ~ celltype, data = veteran,
#'                              dist = "weibull")
#'     celltype_qv2 <- qvcalc(cancer_model2, "celltype")
#'     summary(celltype_qv2)
#' }
#'
#' ## Based on an example from ?itempar
#' if(require(psychotools)) {
#'     data("VerbalAggression", package = "psychotools")
#'     raschmod <- raschmodel(VerbalAggression$resp2)
#'     ip1 <- itempar(raschmod)
#'     qv1 <- qvcalc(ip1)
#'     summary(qv1) }
#'
#' ##  Example of a negative quasi variance
#' ##  Requires the "car" package
#' \dontrun{
#'     library(car)
#'     data(Prestige)
#'     attach(Prestige)
#'     mymodel <- lm(prestige ~ type + education)
#'     library(qvcalc)
#'     type.qvs <- qvcalc(mymodel, "type")
#'     ##  Warning message:
#'     ##  In sqrt(qv) : NaNs produced
#'     summary(type.qvs)
#'     ##  Model call:  lm(formula = prestige ~ type + education)
#'     ##  Factor name:  type
#'     ##          estimate       SE  quasiSE  quasiVar
#'     ##    bc    0.000000 0.000000 2.874361  8.261952
#'     ##    prof  6.142444 4.258961 3.142737  9.876793
#'     ##    wc   -5.458495 2.690667      NaN -1.022262
#'     ##  Worst relative errors in SEs of simple contrasts (%):  0 0
#'     ##  Worst relative errors over *all* contrasts (%):  0 0
#'     plot(type.qvs)
#'     ##  Error in plot.qv(type.qvs) :  No comparison intervals available,
#'     ##  since one of the quasi variances is negative.  See ?qvcalc for more.
#' }
#'
#' @export
qvcalc <- function (object, ...) UseMethod("qvcalc")

#' @export
#' @rdname qvcalc
qvcalc.default <- function(object, factorname = NULL, coef.indices = NULL,
                   labels = NULL, dispersion = NULL,
                   estimates = NULL, modelcall = NULL, ...)
{
    covmat <- object
    if (!is.null(labels)) rownames(covmat) <- colnames(covmat) <- labels
    n <- dim(covmat)[1]
    if (n <= 2) stop(
                    "qvcalc works only for factors with 3 or more levels")
    simple.contrasts <- function(n, levelnames = 1:n){
        result <- list()
        for (i in 1:(n-1)){
            for (j in (i+1):n){
                result[[paste(levelnames[i],
                              levelnames[j],
                              sep = ",")]] <- c(i, j)}}
        result
    }
    qvdesign <- function(n){
        nrows <- choose(n, 2)
        m <- matrix(0, nrows, n)
        indices <- simple.contrasts(n)
        for (i in 1:nrows){
            m[i, indices[[i]][1]] <- 1
            m[i, indices[[i]][2]] <- 1}
        m
    }
    level <- qvdesign(n)
    contrast.variance <- function(contrast, covmat){
        if (!(is.matrix(covmat) &&
              (dim(covmat)[1] == dim(covmat)[2])))
            stop("covmat must be a square matrix")
        n <- dim(covmat)[1]
        if (length(contrast) == n && sum(contrast) == 0)
            ## arbitrary contrast vector
            return(as.vector(contrast %*% covmat %*% contrast))
        if (length(contrast) == 2 && all(contrast %in% 1:n)){
            ## simple contrast specified as an index pair
            i <- contrast[1]
            j <- contrast[2]
            return(covmat[i,i] + covmat[j,j] - 2*covmat[i,j])}
        else stop("invalid contrast")
    }
    simple.contrast.variances <- function(n, covmat){
        if (!is.null(rownames(covmat)))
            levelnames <- rownames(covmat)
        else levelnames <- 1:n
        sapply(simple.contrasts(n, levelnames),
               function(contrast){contrast.variance(contrast, covmat)})
    }
    response <- simple.contrast.variances(n, covmat)
    if (any(response <= 0)) {
        stop("not all contrasts have positive variance")
    } else response <- log(response)
    expLinear <- structure(list(
        family = "expLinear",
        link = "exp",
        linkfun = function(mu) exp(mu),
        linkinv = function(eta) log(eta),
        variance = function(mu) rep(1, length(mu)),
        dev.resids = function(y, mu, wt) wt *
                                         ((y - mu)^2),
        aic = function(y, n, mu, wt, dev) sum(wt) *
                                          (log(dev/sum(wt) * 2 * pi) + 1) + 2,
        mu.eta = function (eta) 1/eta,
        initialize = expression({
            n <- rep(1, nobs)
            mustart <- y}),
        validmu = function(mu) TRUE),
        class = "family")
    model <- glm(response ~ 0 + level, family = expLinear)
    qv <- coef(model)
    NAs <- rep(NA, length(qv))
    if (!is.null(rownames(covmat))) names(qv) <- rownames(covmat)
    frame <- data.frame(estimate = NAs,
                        SE = sqrt(diag(covmat)),
                        quasiSE = sqrt(qv),
                        quasiVar = qv,
                        row.names = names(qv))
    if (!is.null(estimates)) frame$estimate <- estimates
    relerrs <-  sqrt(exp(- residuals(model))) - 1
    ##  The above formula was corrected in v0.8-9; it
    ##  previously said 1 - sqrt(exp(residuals(model)), which is
    ##  not what should be expected for "relative error" here.
    ##  This corrected version agrees with the Biometrika paper.
    ##  Thanks to Shaun Killingbeck for spotting this error in the
    ##  previous version.
    names(relerrs) <- names(response)
    return(structure(list(covmat = covmat,
                          qvframe = frame,
                          dispersion = dispersion,
                          relerrs = relerrs,
                          factorname = factorname,
                          coef.indices = coef.indices,
                          modelcall = modelcall),
                     class="qv"))
}





#' Print with Line Indentation
#'
#' Same as \code{\link{print}}, but adds a specified amount of white space at
#' the start of each printed line
#'
#'
#' @param object any printable object
#' @param indent a non-negative integer, the number of spaces to insert
#' @param \dots other arguments to pass to \code{\link{print}}
#' @return \code{object} is returned invisibly
#' @author David Firth, \email{d.firth@@warwick.ac.uk}
#' @keywords IO
#' @examples
#'
#' indentPrint("this indented by 10 spaces", indent=10)
#'
#' @export
indentPrint <- function(object, indent = 4, ...){
    zz <- ""
    tc <- textConnection("zz", "w", local = TRUE)
    sink(tc)
    try(print(object, ...))
    sink()
    close(tc)
    indent <- paste(rep(" ", indent), sep = "", collapse = "")
    cat(paste(indent, zz, sep = ""), sep = "\n")}

#' @export
print.qv <- function(x, ...){
    print(x$qvframe)
}

#' @export
summary.qv <- function(object, ...)
{
    if (!is.null(object$modelcall))
        cat("Model call: ",
            deparse(object$modelcall), "\n")
    if (!is.null(object$dispersion))
        cat("Dispersion: ", object$dispersion, "\n")
    if (!is.null(object$factorname))
        cat("Factor name: ",object$factorname,"\n")
    indentPrint(object$qvframe,...)
    if (!is.null(object$relerrs)){
        minErrSimple <- round(100*min(object$relerrs), 1)
        maxErrSimple <- round(100*max(object$relerrs), 1)
        errors<-worstErrors(object)
        minErrOverall<-round(100*errors[1], 1)
        maxErrOverall<-round(100*errors[2], 1)
        cat("Worst relative errors in SEs of simple contrasts (%): ",
            minErrSimple, maxErrSimple, "\n")
        cat("Worst relative errors over *all* contrasts (%): ",
            minErrOverall, maxErrOverall, "\n")
    }
}


#' Plot method for objects of class qv
#'
#' Provides visualization of estimated contrasts using intervals based on quasi
#' standard errors.
#'
#' If \code{levelNames} is unspecified, the row names of \code{x$qvframe} will
#' be used.
#'
#' @param x an object of class \code{"qv"}, typically the result of calling
#' \code{\link{qvcalc}}
#' @param intervalWidth the half-width, in quasi standard errors, of the
#' plotted intervals
#' @param ylab as for \code{\link{plot.default}}
#' @param xlab as for \code{\link{plot.default}}
#' @param ylim as for \code{\link{plot.default}}
#' @param main as for \code{\link{plot.default}}
#' @param levelNames labels to be used on the x axis for the levels of the
#' factor whose effect is plotted
#' @param \dots other arguments understood by \code{plot}
#' @return \code{invisible(x)}
#' @author David Firth, \email{d.firth@@warwick.ac.uk}
#' @seealso \code{\link{qvcalc}}
#' @references
#'
#' Easton, D. F, Peto, J. and Babiker, A. G. A. G. (1991) Floating
#' absolute risk: an alternative to relative risk in survival and case-control
#' analysis avoiding an arbitrary reference group.  \emph{Statistics in
#' Medicine} \bold{10}, 1025--1035. \doi{10.1002/sim.4780100703}
#'
#' Firth, D. (2000) Quasi-variances in Xlisp-Stat and on the web.
#' \emph{Journal of Statistical Software} \bold{5.4}, 1--13.
#' \doi{10.18637/jss.v005.i04}
#'
#' Firth, D. (2003) Overcoming the reference category problem in the
#' presentation of statistical models. \emph{Sociological Methodology}
#' \bold{33}, 1--18. \doi{10.1111/j.0081-1750.2003.t01-1-00125.x}
#'
#' Firth, D. and Mezezes, R. X. de (2004) Quasi-variances.  \emph{Biometrika}
#' \bold{91}, 65--80.  \doi{10.1093/biomet/91.1.65}
#'
#' McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
#' London: Chapman and Hall.
#'
#' Menezes, R. X. (1999) More useful standard errors for group and factor
#' effects in generalized linear models.  \emph{D.Phil. Thesis}, Department of
#' Statistics, University of Oxford.
#'
#' @keywords models hplot
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
#' qvs <- qvcalc(shipmodel, "type")
#' summary(qvs, digits = 4)
#' plot(qvs, col = c(rep("red", 4), "blue"))
#' ## if we want to plot in decreasing order (of estimates):
#' est <- qvs$qvframe$estimate
#' qvs2 <- qvs
#' qvs2$qvframe <- qvs$qvframe[order(est, decreasing = TRUE), , drop = FALSE]
#' plot(qvs2)
#'
#' @export
#' 
plot.qv <- function(x,
                    intervalWidth = 2,
                    ylab = "estimate",
                    xlab = "",
                    ylim = NULL,
                    main = "Intervals based on quasi standard errors",
                    levelNames = NULL,
                    ...) {
    frame <- x$qvframe
    if (!is.null(levelNames)) {
        if (nrow(frame) != length(levelNames)) stop(
                                  "levelNames is not a vector of the right length"
                                               )
        row.names(frame) <- levelNames
    }
    if (is.null(frame$quasiSE))
        stop("Cannot plot, because there are no quasi standard errors")
    if (is.na(frame$estimate[1]))
        stop("No parameter estimates to plot")
    if (any(is.nan(frame$quasiSE)))
        stop(paste("No comparison intervals available,\n",
                   "since one of the quasi variances is negative.",
                   "  See ?qvcalc for more.",
                   sep = ""))
    faclevels <- factor(row.names(frame), levels = row.names(frame))
    xvalues <- seq(along = faclevels)
    est <- frame$estimate
    se <- frame$quasiSE
    tops <- est + (intervalWidth * se)
    tails <- est - (intervalWidth * se)
    range <- max(tops) - min(tails)
    if (is.null(ylim)) ylim <- c(min(tails) - range/10, max(tops) + range/10)
    plot(faclevels, frame$estimate, border = "transparent", ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, ...)
    points(frame$estimate, ...)
    segments(xvalues, tails, xvalues, tops, ...)
    invisible(x)
}
