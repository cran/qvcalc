\name{indentPrint}
\alias{indentPrint}
\title{Print with Line Indentation}
\description{
  Same as \code{\link{print}}, but adds a specified amount of white
  space at the start of each printed line
}
\usage{
indentPrint(object, indent=4, ...)
}
\arguments{
  \item{object}{any printable object }
  \item{indent}{a non-negative integer, the number of spaces to insert}
  \item{\dots}{other arguments to pass to \code{\link{print}}}
}
\value{
\code{object} is returned invisibly
}
\author{David Firth, \email{d.firth@warwick.ac.uk}}
\examples{
indentPrint("this indented by 10 spaces", indent=10)
}
\keyword{IO}



