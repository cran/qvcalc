\name{worstErrors}
\alias{worstErrors}
\title{Accuracy of a Quasi-variance Approximation}
\description{
  Computes the worst relative error, among all contrasts,
  for the standard error as derived from a set of quasi-variances.
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
