\name{SIS}
\alias{SIS}
\title{
(Iterative) Sure Independence Screening ((I)SIS) and fitting in Generalized Linear Models and Cox proportional hazards regression models
}
\description{
 This function first implements the iterative sure independence screening with functions \code{\link{GLMvanISISscad}}, \code{\link{GLMvanISISscad}}, \code{\link{COXvanISISscad}}, \code{\link{COXvarISISscad}} 
 for different variants of (I)SIS, then gets the final regression coefficients with functions \code{\link{getfinalSCADcoef}}, \code{\link{INDEPgetfinalSCADcoef}}, 
 \code{\link{getfinalSCADcoefCOX}} for the SCAD regularized loglikelihood for the variables picked by (I)SIS.
 }
\usage{
SIS(data=NULL, model='glm', family=NULL, method='efron', vartype=0, nsis=NULL, 
rank.method='obj', eps0=1e-3, inittype='NoPen', tune.method='BIC', folds=NULL, 
post.tune.method='CV',post.tune.folds=NULL, DOISIS=TRUE, 
ISIStypeCumulative=FALSE, maxloop=5, xtune=NULL, ytune=NULL, detail=FALSE)
}
\arguments{
  \item{data}{
     a list that contains the data.
  }
  \item{model}{
     the model used, the implemented ones are 'glm' and 'cox'.
  }
  \item{family}{
    a description of the error distribution and link function to be used in the generalized linear model.
  }
  \item{method}{
     indicates how to handle observations that have tied (i.e., identical) survival times. The default
    "efron" method is generally preferred to the once-popular "breslow" method.
  }
  \item{vartype}{
   \code{vartype} specifies variant (I)SIS of first type or second type.
   }
  \item{nsis}{
  number of pedictors recuited by (I)SIS.
  }
  \item{rank.method}{
     the criterion for ranking predictor variables in (I)SIS. It can be either \code{obj} or \code{coeff}.
  }
  \item{eps0}{
      an effecitve zero.
  }
  \item{inittype}{
      \code{inittype} specifies the type of initial solution for the one-step SCAD. It can be either \code{NoPen} or \code{L1}.
   }
  \item{tune.method}{
      method for tuning regularization parameter.
  }
  \item{folds}{
    fold information for cross validation.
  }
  \item{post.tune.method}{
    method for tuning regularization parameter in the final step for getting SCAD coefficients.
  }
  \item{post.tune.folds}{
    fold information for cross validation in the final step for getting SCAD coefficients.
  }
  \item{DOISIS}{
   \code{DOISIS} specifies whether to do iterative SIS.
   }
  \item{ISIStypeCumulative}{
   \code{ISIStypeCumulative} specifies whether to allow variable deletion in each step of ISIS. (\code{ISIStypeCumulative}= FALSE allows variable deletion)
   }
   \item{maxloop}{
   maximum number of loops in iterative SIS.
   }
   \item{xtune, ytune}{
   independent tuning dataset.
   }
   \item{detail}{
   indicates whether return detailed information or not. Default is \code{FALSE}.
   }
}
\value{Returns an object with
    \item{SISind}{
    the vector of indices selected by SIS.
    }   
    \item{ISISind}{
     the vector of indices selected by ISIS.
     } 
    \item{SIScoef}{
    a vector of final solution.
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
\code{\link{GLMvanISISscad}}, \code{\link{GLMvanISISscad}}, \code{\link{COXvanISISscad}}, \code{\link{COXvarISISscad}}, \code{\link{getfinalSCADcoef}}, \code{\link{INDEPgetfinalSCADcoef}}, 
 \code{\link{getfinalSCADcoefCOX}}
}
\examples{
set.seed(0)
b <- c(2,2,2,-3*sqrt(2))
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

xtune=matrix(rnorm(n*p, mean=0, sd=1), n, p)
xtune=xtune\%*\%cholmat
feta=xtune[, 1:4]\%*\%b
fprob=exp(feta)/(1+exp(feta))
ytune=rbinom(n, 1, fprob)

binom.result1=SIS(data=list(x=x, y=y), family=binomial(), xtune=xtune, ytune=ytune)
binom.result2=SIS(data=list(x=x, y=y), family=binomial(), xtune=xtune, ytune=ytune, 
vartype=1)
binom.result3=SIS(data=list(x=x, y=y), family=binomial(), xtune=xtune, ytune=ytune, 
vartype=2)



myrates <- exp(x[,1:4]\%*\%b)

ytrue <- rexp(n, rate = myrates) 
cen <- rexp(n, rate = 0.1 )
time <- pmin(ytrue, cen)
status <- as.numeric(ytrue <= cen)

cox.result1=SIS(data=list(x=x,time=time,status=status), model='cox', vartype=0)
cox.result2=SIS(data=list(x=x,time=time,status=status), model='cox', vartype=1)
cox.result3=SIS(data=list(x=x,time=time,status=status), model='cox', vartype=2)


}
\keyword{file}

