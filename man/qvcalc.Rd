\name{qvcalc}
\alias{qvcalc}
\alias{summary.qv}
\alias{plot.qv}
\title{Quasi-variances for Model Coefficients}
\description{
Computes a set of `quasi-variances' (and corresponding `quasi standard errors')
for estimated model coefficients relating to the levels of a categorical
(i.e., factor) explanatory variable.  For details of the method see
Firth (2000) or Firth and Menezes (2002).
}
\usage{
qvcalc(object, factorname=NULL, dispersion = NULL,
      estimates=NULL,  modelcall=NULL)
}
\arguments{
\item{object}{A model (of class lm, glm, etc.), or the covariance (sub)matrix
  for the estimates of interest}
\item{factorname}{If \code{object} is a model, the name of the factor of
  interest}
\item{dispersion}{an optional scalar multiplier for the covariance
  matrix, to cope with overdispersion for example}
\item{estimates}{an optional vector of estimated coefficients (redundant
  if \code{object} is a model)}
\item{modelcall}{optional, the call expression for the model of interest
  (redundant if \code{object} is a model)}
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
\item{modelcall}{if \code{object} is a model, \code{object$call};
  otherwise \code{NULL}}
}
\references{
Firth, D. (2000)  Quasi-variances in Xlisp-Stat and on the web.  
  \emph{Journal of Statistical Software} \bold{5.4}, 1--13.  
  At \url{http://www.jstatsoft.org}

Firth, D. and Mezezes, R. X. de (2002)  Quasi-variances.  Submitted
  for publication.  At \url{http://www.stats.ox.ac.uk/~firth/papers/}.

McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
  London: Chapman and Hall.

Menezes, R. X. (1999)  More useful standard errors for group and factor
  effects in generalized linear models.  \emph{D.Phil. Thesis},
  Department of Statistics, University of Oxford.
}
\author{David Firth, \email{david.firth@nuffield.ox.ac.uk}}

\seealso{\code{\link{worstErrors}}}
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
}
\keyword{models}
\keyword{regression}

