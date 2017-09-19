\name{qvcalc}
\alias{qvcalc}
\alias{qvcalc.default}
\alias{qvcalc.lm}
\alias{qvcalc.itempar}
\alias{summary.qv}
\alias{print.qv}
\title{Quasi Variances for Model Coefficients}
\description{
  Computes a set of quasi variances (and corresponding
  quasi standard errors) for estimated model coefficients
  relating to the levels of a categorical (i.e., factor)
  explanatory variable.  For details of the method see
  Firth (2000), Firth (2003) or Firth and de Menezes (2004).
  Quasi variances generalize and improve the accuracy of
  \dQuote{floating absolute risk} (Easton et al., 1991).
  This device for economical model summary was
  first suggested by Ridout (1989).
}
\usage{
qvcalc(object, ...)

\method{qvcalc}{default}(object, factorname = NULL,  coef.indices = NULL,
      labels = NULL, dispersion = NULL,
      estimates = NULL, modelcall = NULL, ...)

\method{qvcalc}{lm}(object, factorname = NULL, coef.indices = NULL,
       dispersion = NULL, ...)

\method{qvcalc}{itempar}(object, ...)
}
\arguments{
\item{object}{For \code{qvcalc.default}, this is the covariance (sub)matrix
  for the parameters of interest (including any that may have been
  constrained to zero).  For the generic \code{qvcalc}, the
  \code{object} can be any object for which the relevant S3 method has
  been defined.  These currently include many types of regression model
  (via \code{qvcalc.lm}) and objects of class
  \code{\link[psychotools]{itempar}}.}
\item{factorname}{Either \code{NULL}, or a character vector of length 1}
\item{coef.indices}{Either \code{NULL}, or a numeric vector of length at
  least 3}  
\item{labels}{An optional vector of row names for the \code{qvframe}
  component of the result (redundant if \code{object} is a model)}
\item{dispersion}{an optional scalar multiplier for the covariance
  matrix, to cope with overdispersion for example}
\item{estimates}{an optional vector of estimated coefficients (redundant
  if \code{object} is a model, for example)}
\item{modelcall}{optional, the call expression for the model of interest
  (redundant if \code{object} is a model with its own \code{call} component)}
\item{...}{other arguments to pass to \code{qv.default}}
}
\details{

  The \code{qvcalc.default} method is the computational backend for
  all other, class-specific methods.  

  In \code{qvcalc.default}, none of the arguments other than
  \code{object} is used in computing the result.  The remaining
  arguments are simply passed through to the result object as components
  to help with record-keeping etc.

  In \code{qvcalc.lm}, at least one of \code{factorname} or
  \code{coef.indices} must be non-\code{NULL}.  The value of
  \code{coef.indices}, if non-\code{NULL}, determines which rows
  and columns of the model's variance-covariance matrix to use.  If
  \code{coef.indices} contains a zero, then an extra row and column are
  included at the indicated position, to represent the zero variances and
  covariances associated with a reference level.  If \code{coef.indices}
  is \code{NULL}, then \code{factorname} should be the name of a factor
  effect in the model, and is used in order to extract the necessary
  variance-covariance estimates.

  For \code{qvcalc.itempar}, the \code{"itempar"} object must have the
  full variance-covariance matrix in its \code{"vcov"} attribute, and
  must have its \code{"alias"} attribute be \code{TRUE}.  These
  attributes result from use of the default arguments
  \code{vcov = TRUE, alias = TRUE}
  when the \code{\link[psychotools]{itempar}} function is called.
  
  Ordinarily the quasi variances are positive and so their square roots
  (the quasi standard errors) exist and can be used in plots, etc.

  Occasionally one (and only one) of the quasi variances is negative, and
  so the corresponding quasi standard error does not exist (it appears as
  \code{NaN}).  This is fairly rare in applications, and
  when it occurs it is because the factor of interest is strongly
  correlated with one or more other predictors in the model.  It is not
  an indication that quasi variances are inaccurate.  An example is shown
  below using
  data from the \code{car} package: the quasi variance
  approximation
  is exact (since \code{type} has
  only 3 levels), and there is a negative quasi variance.  The quasi
  variances remain perfectly valid (they can be used to obtain
  inference on any contrast), but it makes no sense to plot
  `comparison intervals' in the usual way since one of the quasi standard
  errors is not a real number.
}
\value{
A list of class \code{qv}, with components
\item{covmat}{the full variance-covariance matrix for the estimated 
  coefficients corresponding to the factor of interest}
\item{qvframe}{a data frame with variables
  \code{estimate}, \code{SE}, \code{quasiSE} and \code{quasiVar},
  the last two being a quasi standard error and quasi-variance
  for each level of the factor of interest}
\item{relerrs}{relative errors for approximating the standard errors of all 
  simple contrasts}
\item{factorname}{the factor name if given}
\item{coef.indices}{the coefficient indices if given}
\item{modelcall}{if \code{object} is a model, \code{object$call};
  otherwise \code{NULL}}
}
\references{
Easton, D. F, Peto, J. and Babiker, A. G. A. G. (1991)  
Floating absolute risk: an alternative to relative risk in survival 
and case-control analysis avoiding an arbitrary reference group.  
\emph{Statistics in Medicine} \bold{10}, 1025--1035.

Firth, D. (2000)  Quasi-variances in Xlisp-Stat and on the web.  
  \emph{Journal of Statistical Software} \bold{5.4}, 1--13.  
  At \url{http://www.jstatsoft.org}

Firth, D. (2003) Overcoming the reference category problem in the 
presentation of statistical models. \emph{Sociological Methodology} 
\bold{33}, 1--18.

Firth, D. and de Mezezes, R. X. (2004)  Quasi-variances.  
  \emph{Biometrika} \bold{91}, 65--80.  

McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
  London: Chapman and Hall.

Menezes, R. X. de (1999)  More useful standard errors for group and factor
  effects in generalized linear models.  \emph{D.Phil. Thesis},
  Department of Statistics, University of Oxford.

Ridout, M.S. (1989). Summarizing the results of fitting generalized
  linear models to data from designed experiments. In: \emph{Statistical
  Modelling: Proceedings of GLIM89 and the 4th International
  Workshop on Statistical Modelling held in Trento, Italy, July 17--21,
  1989} (A. Decarli et al., eds.), pp 262--269. New York: Springer.
  
}
\author{David Firth, \email{d.firth@warwick.ac.uk}}

\seealso{
  \code{\link{worstErrors}}, \code{\link{plot.qv}}.  
}

\examples{
##  Overdispersed Poisson loglinear model for ship damage data
##  from McCullagh and Nelder (1989), Sec 6.3.2 
if (require(MASS)) {
    data(ships)
    ships$year <- as.factor(ships$year)
    ships$period <- as.factor(ships$period)
    shipmodel <- glm(formula = incidents ~ type + year + period,
                     family = quasipoisson, 
                     data = ships,
                     subset = (service > 0),
                     offset = log(service))
    shiptype.qvs <- qvcalc(shipmodel, "type")
    summary(shiptype.qvs, digits = 4)
    plot(shiptype.qvs, xlab = "ship type")
## An equivalent result by using the coef.indices argument instead:
    shiptype.qv2 <- qvcalc(shipmodel, coef.indices = c(0, 2:5))
}

## Based on an example from ?itempar
if(require(psychotools)) {
    data("VerbalAggression", package = "psychotools")
    raschmod <- raschmodel(VerbalAggression$resp2)
    ip1 <- itempar(raschmod)
    qv1 <- qvcalc(ip1)
    summary(qv1) }

##  Example of a negative quasi variance
##  Requires the "car" package
\dontrun{
    library(car)
    data(Prestige)
    attach(Prestige)
    mymodel <- lm(prestige ~ type + education)
    library(qvcalc)
    type.qvs <- qvcalc(mymodel, "type")
    ##  Warning message: 
    ##  In sqrt(qv) : NaNs produced
    summary(type.qvs)
    ##  Model call:  lm(formula = prestige ~ type + education) 
    ##  Factor name:  type 
    ##          estimate       SE  quasiSE  quasiVar
    ##    bc    0.000000 0.000000 2.874361  8.261952
    ##    prof  6.142444 4.258961 3.142737  9.876793
    ##    wc   -5.458495 2.690667      NaN -1.022262
    ##  Worst relative errors in SEs of simple contrasts (\%):  0 0 
    ##  Worst relative errors over *all* contrasts (\%):  0 0
    plot(type.qvs)
    ##  Error in plot.qv(type.qvs) :  No comparison intervals available,
    ##  since one of the quasi variances is negative.  See ?qvcalc for more.
} 
}
\keyword{models}
\keyword{regression}
