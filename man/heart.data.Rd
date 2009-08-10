\name{heart.data}
\alias{heart.data}
\docType{data}
\title{Dataset for SIS}
\description{
  \emph{South African Heart Disease dataset} used to test \code{SIS}
  algorithm
}
\usage{
data(heart.data)
}
\format{
  A dataset with 462 observations on 9 variables and a binary
  response.
  \describe{
    \item{x}{x contains 9 columns of the following variables:
      sbp (systolic blood pressure); tobacco (cumulative tobacco); ldl
      (low density lipoprotein cholesterol); adiposity; famhist (family
      history of heart disease); typea (type-A behavior); obesity;
      alcohol (current alcohol consumption); age (age at onset)}
    \item{y}{response, coronary heart disease}
    }
}
\examples{
data(heart.data)
attach(heart.data)
n=length(y)
p=dim(x)[2]
set.seed(0)
trainset=sample(1:n,floor(n*2/3))
testset=setdiff(1:n,trainset)

traindata=list(x=x[trainset,],y=y[trainset])

fit <- SIS(traindata, family=binomial(), xtune=x[testset,], ytune=y[testset])

detach(heart.data)
}
\references{
  Hastie, T., Tibshirani, R., and Friedman, J. (2001)
  \emph{Elements of Statistical Learning; Data Mining, Inference, and
    Prediction} Springer-Verlag, New York.
}
\keyword{datasets}
