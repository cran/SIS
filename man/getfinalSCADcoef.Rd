\name{getfinalSCADcoef}
\alias{getfinalSCADcoef}
\alias{INDEPgetfinalSCADcoef}
\title{
  SCAD regularized loglikelihood for generalized linear models
  
}
\description{
 This function gets the final regression coefficients for the SCAD regularized loglikelihood for generalized linear models after applying (I)SIS.
  }
\usage{
getfinalSCADcoef(x, y, pickind, folds=NULL, eps0=1e-3, family=binomial(), 
tune.method="AIC", inittype='NoPen', detailed=FALSE)

INDEPgetfinalSCADcoef(x, y, pickind,  folds=NULL, xtune, ytune, 
family=binomial(), inittype='NoPen', eps0=1e-3, detailed=FALSE) 

}
\arguments{
  \item{x}{
     an (n * p)  matrix of features.
  }
  \item{y}{
     an (n) vector of response.
  }
  \item{pickind}{
    predictor indice selected by (I)SIS.
  }
  \item{folds}{
    fold information for cross validation.
  }
  \item{family}{
    a description of the error distribution and link function to be used in the model.
    }
  \item{tune.method}{
      method for tuning regularization parameter.
    }
  \item{inittype}{
      \code{inittype} specifies the type of initial solution for the one-step SCAD. It can be either \code{NoPen} or \code{L1}.
   }
   \item{xtune, ytune}{
   independent tuning dataset.
   }
  \item{eps0}{
      an effective zero.
  }
  \item{detailed}{
    indicates whether detailed information will be provided. Default is \code{FALSE}.
  }
}
\value{An initial solution vector \code{wt.initsoln} and final solution (p+1) vector \code{SCADcoef}.
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
\code{\link{scadglm}}, \code{\link{fullscadglm}}
}
\examples{
set.seed(0)
b <- c(4,4,4,-6*sqrt(2))
n=150
p=200
truerho=0.5
corrmat=diag(rep(1-truerho, p))+matrix(truerho, p, p)
corrmat[,4]=sqrt(truerho)
corrmat[4, ]=sqrt(truerho)
corrmat[4,4]=1
cholmat=chol(corrmat)
x=matrix(rnorm(n*p, mean=0, sd=1), n, p)
x=x\%*\%cholmat
feta=x[, 1:4]\%*\%b
fprob=exp(feta)/(1+exp(feta))
y=rbinom(n, 1, fprob)

nsis=floor(n/log(n)/4)
xtune=matrix(rnorm(n*p, mean=0, sd=1), n, p)
xtune=xtune\%*\%cholmat
feta=xtune[, 1:4]\%*\%b
fprob=exp(feta)/(1+exp(feta))
ytune=rbinom(n, 1, fprob)

ISIScoef = INDEPgetfinalSCADcoef(x = x, y = y, 
                pickind = 1:4, xtune = xtune, ytune = ytune, 
                family = binomial())
}

