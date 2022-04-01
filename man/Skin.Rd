\name{Skin}
\alias{Skin}
\docType{data}
\title{
Skin segmentation
}
\description{
The skin segmentation data
}
\usage{data("Skin")}
\format{
  A data frame with 245057 observations on the following 4 variables.
  \describe{
    \item{\code{B}}{a numeric vector}
    \item{\code{G}}{a numeric vector}
    \item{\code{R}}{a numeric vector}
    \item{\code{C}}{a numeric vector}
  }
}
\details{
The skin segmentation data is related to skin texture in face image. The total number of samples is 245057, and the feature number is 3.
}
\source{
The skin segmentation data set is from the UCI database. 
}
\references{
Rajen B. Bhatt, Gaurav Sharma, Abhinav Dhall, Santanu Chaudhury, Efficient skin region segmentation using low complexity fuzzy decision tree model, IEEE-INDICON 2009, Dec 16-18, Ahmedabad, India, pp. 1-4.
}
\examples{
data(Skin)
## maybe str(Skin) ; plot(Skin) ...
}
\keyword{datasets}
