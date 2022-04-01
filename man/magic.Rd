\name{magic}
\alias{magic}
\docType{data}
\title{
Magic
}
\description{
The magic data
}
\usage{data("magic")}
\format{
  A data frame with 19020 observations on the following 11 variables.
  \describe{
    \item{\code{fLength}}{a numeric vector}
    \item{\code{fWidth}}{a numeric vector}
    \item{\code{fSize}}{a numeric vector}
    \item{\code{fConc}}{a numeric vector}
    \item{\code{fConc1}}{a numeric vector}
    \item{\code{fAsym}}{a numeric vector}
    \item{\code{fM3Long}}{a numeric vector}
    \item{\code{fM3Trans}}{a numeric vector}
    \item{\code{fAlpha}}{a numeric vector}
    \item{\code{fDist}}{a numeric vector}
    \item{\code{class}}{a character vector}
  }
}
\details{
The magic data set is given by MAGIC project, and described by 11 features.
}
\source{
The magic data set is from the UCI database. 
}
\references{
J. Dvorak, P. Savicky.  Softening Splits in Decision Trees Using Simulated Annealing.  Proceedings of ICANNGA 2007, Warsaw, Part I, LNCS 4431, pp. 721-729.
}
\examples{
data(magic)
## maybe str(magic) ; plot(magic) ...
}
\keyword{datasets}
