\name{wtlassoglm}
\alias{wtlassoglm}
\title{
  Weighted L1 regularized loglikelihood for generalized linear models
}
\description{
  This functions solves weighted L1 regularized loglikelihood for generalized linear models.
}
\usage{
wtlassoglm(x, y, lassoweight=NULL, initsoln=NULL,   family = binomial(), 
weight = NULL, offset = NULL, lambda2=0, function.precision=1e-8)
}
\arguments{
  \item{x}{
     an (n * p)  matrix of features. 
  }
  \item{y}{
     an (n) vector of response.
  }
  \item{lassoweight}{
    a (p) vector of weights specifying the weighted L1 penalty.
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
      this can be used to specify an \code{a priori} known component to be included 
      in the linear predictor during fitting.
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
  \item{family}{
    a description of the error distribution and link function to be used in the model.
    }
  \item{weight}{
      an optional (n) vector of weights to be used in the fitting process.
    }
  \item{offset}{
      this can be used to specify an \code{a priori} known component to be included in 
      the linear predictor during fitting.
   }
  \item{lassoweight}{
    a (p) vector of weights specifying the weighted L1 penalty.
  }
  \item{initsoln}{
    a (p+1) vector of initial solution. 
  }
  \item{w}{
    a (p+1) vector of weight L1 solution.
    }
}
\details{
This function solves weighted L1 regularized loglikelihood for generalized linear models.
It is based on the source code of R package \code{glmpath}.
}

\author{Jianqing Fan, Yang Feng, Richard Samworth, and Yichao Wu}
\seealso{
\code{\link{scadglm}}, \code{\link{fullscadglm}}
}
\examples{
set.seed(0)
b <- c(2,2,2,-3*sqrt(2))
n=400
p=30
truerho=0.5
x=matrix(rnorm(n*p, mean=0, sd=1), n, p)
feta=x[, 1:4]\%*\%b
fprob=exp(feta)/(1+exp(feta))
y=rbinom(n, 1, fprob)
lassoweight<-rep(0.6,30)
wtlassoglm(x,y,lassoweight)$w
coef(glm(y~x,family=binomial()))
}

\keyword{file}
