\name{wtlassocox}
\alias{wtlassocox}
\title{
  Weighted L1 regularized loglikelihood for Cox proportional hazards regression model
}
\description{
  This functions solves weighted L1 regularized loglikelihood for Cox proportional hazards regression model.
}
\usage{
wtlassocox(x, time, status, method="efron", lassoweight=NULL, initsoln=NULL,  
weight = NULL, lambda2=0, function.precision=1e-8)
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
    
  \item{lassoweight}{
    an optional (p) vector of weights specifying the weighted L1 penalty.
  }
  \item{initsoln}{
    an optional (p) vector of initial solution.
  }
  \item{weight}{
    an optional (n) vector of weights to be used in the fitting process.
  }
  \item{lambda2}{
    regularization parameter for the L2 norm of the
    coefficients. Default is \code{0.}
  }
  \item{function.precision}{
    \code{function.precision} parameter used in the internal
    solver. Default is \code{1e-8.}
  }
}
\value{
  An object is returned with
  \item{lambda2}{
    \eqn{\lambda_2} used.
  }
  \item{xnames}{
    column names of x.
  }
  \item{weight}{
      an optional (n) vector of weights to be used in the fitting process.
    }
  \item{lassoweight}{
    a (p) vector of weights specifying the weighted L1 penalty.
  }
  \item{initsoln}{
    a (p) vector of initial solution.
  }
  \item{w}{
    a (p) vector of weight L1 solution.
    }
}
\details{
This function solves weighted L1 regularized loglikelihood for Cox proportional hazards regression model.
It is based on the source code of R package \code{glmpath}.
}

\author{Jianqing Fan, Yang Feng, Richard Samworth, and Yichao Wu}
\seealso{
\code{\link{scadglm}}, \code{\link{fullscadglm}}
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



weights <- rep(0.01*n, 100)
w1 <- wtlassocox(x[,1:100],time,status,lassoweight=weights)$w
w2 <- coxph(Surv(time,status)~x[,1:100])$coef
}

