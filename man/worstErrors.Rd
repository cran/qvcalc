\name{worstErrors}
\alias{worstErrors}
\title{Accuracy of a Quasi-variance Approximation}
\description{
  Computes the worst relative error, among all contrasts,
  for the standard error as derived from a set of quasi-variances.
  For details of the method see Menezes (1999) or Firth and
  Menezes (2002).
}
\usage{
worstErrors(qv.object)
}
\arguments{
  \item{qv.object}{An object of class \code{qv}}
}
\value{
  A numeric vector of length 2, the worst negative relative
  error and the worst positive relative error.
}
\references{
Firth, D. and Mezezes, R. X. de (2002)  Quasi-variances.  Submitted
  for publication.  At \url{http://www.stats.ox.ac.uk/~firth/papers/}.

McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
  London: Chapman and Hall.

Menezes, R. X. (1999)  More useful standard errors for group and factor
  effects in generalized linear models.  \emph{D.Phil. Thesis},
  Department of Statistics, University of Oxford.
}
\author{David Firth, \email{david.firth@nuffield.ox.ac.uk}}
\seealso{\code{\link{qvcalc}}}

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
worstErrors(shiptype.qvs)
}

\keyword{regression}
\keyword{models}
