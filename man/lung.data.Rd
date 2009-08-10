\name{lung.data}
\alias{lung.data}
\docType{data}
\title{Dataset for SIS}
\description{
  \emph{Lung cancer dataset} used to test \code{SIS} algorithm
}
\usage{
data(lung.data)
}
\format{
  A dataset consisting of 137 observations with their survival time,
  censor status as well as 6 features. 
  \describe{
    \item{x}{x contains 6 columns of the following variables: trt
      (1=standard treatment, and 2=test); celltype (1=squamous,
      2=smallcell, 3=adeno, and 4=large); karno (Karnofsky performance
      score); diagtime (months from diagnosis to randomization); age (in
      years); prior (prior therapy 0=no, and 1=yes)}
    \item{time}{survival time}
    \item{status}{censor status}
    }
}
\examples{

data(lung.data)
attach(lung.data)
fit <- SIS(lung.data, model='cox')
detach(lung.data)
}
\references{
  Kalbfleisch, J. and Prentice, R. (2002) \emph{The Statistical Analysis
  of Failure Time Data} J. Wiley, Hoboken, N.J.
}
\keyword{datasets}
