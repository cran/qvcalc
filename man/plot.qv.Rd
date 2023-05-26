\name{plot.qv}
\alias{plot.qv}
\title{ Plot method for objects of class qv }
\description{
  Provides visualization of estimated contrasts using intervals based on
  quasi standard errors.
}
\usage{
\method{plot}{qv}(x, intervalWidth = 2, ylab = "estimate",
    xlab = "", ylim = NULL,
    main = "Intervals based on quasi standard errors",
    levelNames = NULL, ...)
}
\arguments{
  \item{x}{ an object of class \code{"qv"}, typically the result of
  calling \code{\link{qvcalc}}}
  \item{intervalWidth}{ the half-width, in quasi standard errors, of the
  plotted intervals}
  \item{ylab}{ as for \code{\link{plot.default}}}
  \item{xlab}{ as for \code{\link{plot.default}}}
  \item{ylim}{ as for \code{\link{plot.default}}}
  \item{main}{ as for \code{\link{plot.default}}}
  \item{levelNames}{ labels to be used on the x axis for the levels of
  the factor whose effect is plotted }
  \item{\dots}{ other arguments understood by \code{plot} }
}
\details{
  If \code{levelNames} is unspecified, the row names of \code{x$qvframe}
  will be used.

}
\value{
  \code{invisible(x)}
}
\references{
Easton, D. F, Peto, J. and Babiker, A. G. A. G. (1991)  
Floating absolute risk: an alternative to relative risk in survival 
and case-control analysis avoiding an arbitrary reference group.  
\emph{Statistics in Medicine} \bold{10}, 1025--1035.

Firth, D. (2000)  Quasi-variances in Xlisp-Stat and on the web.  
  \emph{Journal of Statistical Software} \bold{5.4}, 1--13.  
  \doi{10.18637/jss.v005.i04}

Firth, D. (2003) Overcoming the reference category problem in the 
presentation of statistical models. \emph{Sociological Methodology} 
\bold{33}, 1--18. \doi{10.1111/j.0081-1750.2003.t01-1-00125.x}

Firth, D. and Mezezes, R. X. de (2004)  Quasi-variances.  
  \emph{Biometrika} \bold{91}, 65--80.  \doi{10.1093/biomet/91.1.65}

McCullagh, P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}.
  London: Chapman and Hall.

Menezes, R. X. (1999)  More useful standard errors for group and factor
  effects in generalized linear models.  \emph{D.Phil. Thesis},
  Department of Statistics, University of Oxford.
}
\author{David Firth, \email{d.firth@warwick.ac.uk}}
\seealso{ \code{\link{qvcalc}} }
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
qvs <- qvcalc(shipmodel, "type")
summary(qvs, digits = 4)
plot(qvs, col = c(rep("red", 4), "blue"))
## if we want to plot in decreasing order (of estimates):
est <- qvs$qvframe$estimate
qvs2 <- qvs
qvs2$qvframe <- qvs$qvframe[order(est, decreasing = TRUE), , drop = FALSE]
plot(qvs2)
}
\keyword{ models }
\keyword{ hplot }
