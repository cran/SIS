\name{COXvanISISscad}
\alias{COXvanISISscad}
\alias{COXvarISISscad}
\title{
(Iterative) Sure Independence Screening ((I)SIS) in the Cox proportional hazards regression model
}
\description{
 These functions implement the iterative sure independence screening with \code{\link{COXvanISISscad}} for vanilla ISIS and  \code{\link{COXvarISISscad}} for variant ISIS in Cox proportional hazards regression model.
  }
\usage{
COXvanISISscad(x, time, status,  method = "efron", nsis=NULL, folds=NULL, 
rank.method="obj", eps0=1e-3, inittype='NoPen', tune.method="AIC", 
ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5)

COXvarISISscad(x, time, status,  method = "efron", nsis=NULL, folds=NULL, 
rank.method="obj", eps0=1e-3, inittype='NoPen', tune.method="AIC", 
vartype="First",  ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5)
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
  \item{nsis}{
  number of pedictors recuited by (I)SIS.
  }
  \item{folds}{
    fold information for cross validation.
  }
  \item{rank.method}{
     the criterion for ranking predictor variables in (I)SIS. It can be either \code{obj} or \code{coeff}.
     }
  \item{tune.method}{
      method for tuning regularization parameter.
    }
  \item{inittype}{
      \code{inittype} specifies the type of initial solution for the one-step SCAD. It can be either \code{NoPen} or \code{L1}.
   }
   \item{vartype}{
   \code{vartype} specifies variant (I)SIS of first type or second type.
   }
   \item{ISIStypeCumulative}{
   \code{ISIStypeCumulative} specifies whether to allow variable deletion in each step of ISIS. (\code{ISIStypeCumulative}= FALSE allows variable deletion)
   }
   \item{DOISIS}{
   \code{DOISIS} specifies whether to do iterative SIS.
   }
   \item{maxloop}{
   maximum number of loops in iterative SIS.
   }
  \item{eps0}{
      an effective zero.
  }
}
\value{Returns an object with
\item{initRANKorder}{
    initial predictor ranking order for vanilla SIS.
    }
\item{detail.pickind, detail.ISISind}{
   details of each loop of ISIS.
   }
\item{SISind}{
 the vector of indices selected by SIS.
  }
\item{ISISind}{
 the vector of indices selected by ISIS.
 }
 \item{normal.exit}{
   indicator of normal exit.
   }
 \item{nsis}{
 number of pedictors recuited by (I)SIS.
 }
\item{initRANKorder1, initRANKorder2}{
   initial predictor ranking order for variant SIS.
   }
    


}

\references{
Jianqing Fan and Jinchi Lv (2008) Sure independence screening for ultra-high dimensional feature space
(with discussion) \emph{Journal of Royal Statistical Society B}, \bold{36}, 849-911.

Jianqing Fan, Richard Samworth, and Yichao Wu (2009) Ultrahigh dimensional variable selection: beyond the linear
model \emph{Journal of Machine Learning Research}, to appear.

Jianqing Fan and Rui Song (2009) Sure Independence Screening in Generalized Linear Models with NP-Dimensionality, technical report.

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

cox.van.sis <- COXvanISISscad(x, time, status)

cox.var.sis <- COXvarISISscad(x, time, status, vartype='Second')

#####compare the result
cox.van.sis$SIS

cox.van.sis$ISIS

cox.var.sis$SIS

cox.var.sis$ISIS


}

