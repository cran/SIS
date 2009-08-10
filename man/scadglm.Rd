\name{scadglm}
\alias{scadglm}
\alias{fullscadglm}
\title{
  SCAD regularized loglikelihood for generalized linear models
}
\description{
These functions solve SCAD regularized loglikelihood for generalized linear models; \code{scadcox} does the
    one-step SCAD while \code{fullscadcox} solves the SCAD in a fully
    iterative method.
}
\usage{
scadglm(x, y, wt.initsoln=NULL, lambda, initsoln=NULL, 
family = binomial(), weight = NULL, offset = NULL, 
function.precision=1e-8, nopenalty.subset=NULL)
       
fullscadglm(x, y, lambda, initsoln=NULL, family = binomial(),
   weight = NULL, offset = NULL, function.precision=1e-8, 
   nopenalty.subset=NULL, eps0=1e-6)

}
\arguments{
  \item{x}{
     an (n * p)  matrix of features.
  }
  \item{y}{
     an (n) vector of response.
  }
  \item{wt.initsoln}{
    a (p+1) vector of initial solution for one-step SCAD.
  }
  \item{lambda}{
    regularization parameter for the SCAD.
  }
  \item{initsoln}{
    a (p+1) vector of initial solution.
  }
  \item{family}{
    a description of the error distribution and link function to be used in the model.
    }
  \item{weight}{
      an optional (n) vector of weights to be used in the fitting process.
    }
  \item{offset}{
      this can be used to specify an \code{a priori} known component to be included in the linear predictor during fitting.
   }
  \item{function.precision}{
    \code{function.precision} parameter used in the internal
    solver. Default is \code{1e-8.}
  }
  \item{nopenalty.subset}{
    a set of indices for the predictors that are not subject to the L1
    penalty.
  }
  \item{eps0}{
      an effective zero.
  }
}
\value{  They return a (p+1) vector of estimated coefficients.
}

\references{
 Jianqing Fan and Runze Li (2001) Variable selection via nonconcave
 penalized likelihood and its oracle properties. \emph{Journal of American
  Statistical Association}, \bold{96}, 1348-1360.
  
 Hui Zou and Runze Li (2008) One-step Sparse Estimates in Nonconcave
 Penalized Likelihood Models (with discussion). \emph{The Annals of Statistics}, \bold{36}, 1509-1533
  
  }
\author{Jianqing Fan, Yang Feng, Richard Samworth, and Yichao Wu}
\seealso{
\code{\link{wtlassoglm}}
}
\examples{
set.seed(0)
b <- c(1,1,1,-3*sqrt(2)/2)
n=400
p=30
truerho=0.5
x=matrix(rnorm(n*p, mean=0, sd=1), n, p)
feta=x[, 1:4]\%*\%b
fprob=exp(feta)/(1+exp(feta))
y=rbinom(n, 1, fprob)
scadglm(x,y,lambda=0.0015)
coef(glm(y~x,family=binomial()))
}
