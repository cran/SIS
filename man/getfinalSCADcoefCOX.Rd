\name{getfinalSCADcoefCOX}
\alias{getfinalSCADcoefCOX}
\title{
  SCAD regularized loglikelihood for Cox proportional hazards regression models

}
\description{
 This function gets the final regression coefficients for the SCAD regularized loglikelihood for Cox proportional hazards regression models after applying (I)SIS
  }
\usage{
getfinalSCADcoefCOX(x, time, status, method = "efron", pickind, 
folds = NULL, eps0 = 1e-3, tune.method = "AIC", inittype = "NoPen", 
detailed = FALSE)
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
  \item{pickind}{
    predictor indice selected by (I)SIS.
  }
  \item{folds}{
    fold information for cross validation.
  }
  \item{eps0}{
      an effecitve zero.
  }
  \item{tune.method}{
      method for tuning regularization parameter.
    }
  \item{inittype}{
      \code{inittype} specifies the type of initial solution for the one-step SCAD. It can be either \code{NoPen} or \code{L1}.
   }
  \item{detailed}{
    indicates whether detailed information will be provided. Default is \code{FALSE}.
  }
}
\value{
An initial solution vector \code{wt.initsoln} and final solution (p) vector \code{SCADcoef}.
}

\references{
Jianqing Fan and Runze Li (2002) Variable Selection for Cox's Proportional 
 Hazards Model and Frailty Model. \emph{The Annals of Statistics}, \bold{30}, 74-99.

Hui Zou and Runze Li (2008) One-step Sparse Estimates in Nonconcave
 Penalized Likelihood Models (with discussion). \emph{The Annals of Statistics}, \bold{36}, 1509-1533
}
\author{Jianqing Fan, Yang Feng, Richard Samworth, and Yichao Wu}
\seealso{
\code{\link{scadcox}}, \code{\link{fullscadcox}}
}
\examples{
set.seed(0)
n=150
p=200
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

SIScoef <- getfinalSCADcoefCOX(x = x, time = time, status = status, 
            pickind = 1:5)
}

\keyword{file}
