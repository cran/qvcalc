\name{qvcalc}
\alias{qvcalc}
\alias{summary.qv}
\alias{print.qv}
\title{Quasi Variances for Model Coefficients}
\description{
  Computes a set of quasi variances (and corresponding
  quasi standard errors)
for estimated model coefficients relating to the levels of a categorical
(i.e., factor) explanatory variable.  For details of the method see
Firth (2000), Firth (2003) or Firth and Menezes (2004).  Quasi variances 
generalize
and improve the accuracy of \dQuote{floating absolute risk} 
(Easton et al., 1991). 
}
\usage{
qvcalc(object, factorname=NULL,  coef.indices = NULL,
      labels = NULL, dispersion = NULL,
      estimates=NULL,  modelcall=NULL)
}
\arguments{
\item{object}{A model (of class lm, glm, etc.), or the covariance (sub)matrix
  for the estimates of interest, or an object of class \code{Btabilities}}
\item{factorname}{Either \code{NULL}, or a character vector of length 1}
\item{coef.indices}{Either \code{NULL}, or a numeric vector of length at
  least 3}  
\item{labels}{An optional vector of row names for the \code{qvframe}
  component of the result (redundant if \code{object} is a model)}
\item{dispersion}{an optional scalar multiplier for the covariance
  matrix, to cope with overdispersion for example}
\item{estimates}{an optional vector of estimated coefficients (redundant
  if \code{object} is a model)}
\item{modelcall}{optional, the call expression for the model of interest
  (redundant if \code{object} is a model)}
}
\details{
If \code{object} is a model, then at least one of \code{factorname} or
\code{coef.indices} must be non-\code{NULL}.  The value of
\code{coef.indices}, if non-\code{NULL}, determines which rows
and columns of the model's variance-covariance matrix to use.  If
\code{coef.indices} contains a zero, an extra row and column are
included at the indicated position, to represent the zero variances and
covariances associated with a reference level.  If \code{coef.indices}
is \code{NULL}, then \code{factorname} should be the name of a factor
effect in the model, and is used in order to extract the necessary
variance-covariance estinmates.

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

Firth, D. and Mezezes, R. X. de (2004)  Quasi-variances.  
  \emph{Biometrika} \bold{91}, 65--80.  

McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
  London: Chapman and Hall.

Menezes, R. X. (1999)  More useful standard errors for group and factor
  effects in generalized linear models.  \emph{D.Phil. Thesis},
  Department of Statistics, University of Oxford.
}
\author{David Firth, \email{d.firth@warwick.ac.uk}}

\seealso{\code{\link{worstErrors}}, \code{\link{plot.qv}}}
\examples{
##  Overdispersed Poisson loglinear model for ship damage data
##  from McCullagh and Nelder (1989), Sec 6.3.2 
library(MASS)
data(ships)
ships$year <- as.factor(ships$year)
ships$period <- as.factor(ships$period)
shipmodel <- glm(formula = incidents ~ type + year + period,
    family = quasipoisson, 
    data = ships, subset = (service > 0), offset = log(service))
shiptype.qvs <- qvcalc(shipmodel, "type")
summary(shiptype.qvs, digits=4)
plot(shiptype.qvs)

##  Quasi-variance summary for "ability" estimates in a Bradley-Terry model
##  Requires the "BradleyTerry2" package
\dontrun{
library(BradleyTerry2)
example(baseball)
plot(baseball.qv <- qvcalc(BTabilities(baseballModel2)))
}


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
##  NaNs produced in: sqrt(qv) 
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

