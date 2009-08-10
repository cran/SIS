\name{scadcox}
\alias{scadcox}
\alias{fullscadcox}
\title{
  SCAD regularized loglikelihood for Cox proportional hazards regression models
}
\description{
These functions solve SCAD regularized loglikelihood for Cox
    proportional hazards regression models; \code{scadcox} does the
    one-step SCAD while \code{fullscadcox} solves the SCAD in a fully
    iterative method.
}
\usage{
scadcox(x, time, status,  method="efron", wt.initsoln=NULL, lambda, 
    initsoln=NULL, weight = NULL, function.precision=1e-8, 
    nopenalty.subset=NULL)
    
fullscadcox(x, time, status,  method="efron", lambda, 
    initsoln=NULL, weight = NULL, function.precision=1e-8, 
    nopenalty.subset=NULL, eps0=1e-6)
}
\arguments{
  \item{x}{
    an (n * p)  matrix of features.
  }
  \item{time}{
    an (n) vector of the follow up time for right censored data.
  }
  \item{status}{
    an (n) vector of the status indicator, normally 0=alive, 1=dead.
  }
  \item{method}{
     indicates how to handle observations that have tied (i.e., identical) survival times. The default
    "efron" method is generally preferred to the once-popular "breslow" method.
  }
  \item{wt.initsoln}{
    a (p) vector of initial solution for one-step SCAD.
  }
  \item{lambda}{
    the regularization parameter.
  }
  \item{initsoln}{
    a (p) vector of initial solution.
  }
  \item{weight}{
      an optional (n) vector of weights to be used in the fitting process.
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
    an effecitve zero.
  }
}
\value{They return a (p) vector of estimated coefficients.
}

\references{
Jianqing Fan and Runze Li (2002) Variable Selection for Cox's Proportional 
 Hazards Model and Frailty Model. \emph{The Annals of Statistics}, \bold{30}, 74-99.

 Hui Zou and Runze Li (2008) One-step Sparse Estimates in Nonconcave
 Penalized Likelihood Models (with discussion). \emph{The Annals of Statistics}, \bold{36}, 1509-1533

  }
\author{Jianqing Fan, Yang Feng, Richard Samworth, and Yichao Wu}
\seealso{
\code{\link{wtlassocox}}
}
\examples{
set.seed(0)
n=400
p=1000
truerho=0.5
beta <- c(4,4,4,-6*sqrt(2),4/3, rep(0,p-5))

corrmat=diag(rep(1-truerho, p))+matrix(truerho, p, p)
corrmat[,4]=sqrt(truerho)
corrmat[4, ]=sqrt(truerho)
corrmat[4,4]=1
corrmat[,5]=0
corrmat[5,]=0
corrmat[5,5]=1
cholmat=chol(corrmat)

x=matrix(rnorm(p*n, mean=0, sd=1), n, p)
x=x\%*\%cholmat

myrates <- exp(x\%*\%beta)

ytrue <- rexp(n, rate = myrates) 
cen <- rexp(n, rate = 0.1 )
time <- pmin(ytrue, cen)
status <- as.numeric(ytrue <= cen)

w1 <- scadcox(x[,1:100],time,status,lambda=0.003)
w2 <- coxph(Surv(time,status)~x[,1:100])$coef


}

