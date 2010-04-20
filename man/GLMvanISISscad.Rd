\name{GLMvanISISscad}
\alias{GLMvanISISscad}
\alias{GLMvarISISscad}
\title{
(Iterative) Sure Independence Screening ((I)SIS) in Generalized Linear Models
}
\description{
 These functions implement the iterative sure independence screening with \code{\link{GLMvanISISscad}} for vanilla ISIS and  \code{\link{GLMvarISISscad}} for variant ISIS in Generalized Linear Models.
  }
\usage{
GLMvanISISscad(x, y, nsis=NULL, family=binomial(), folds=folds, 
rank.method="obj", eps0=1e-5, inittype='NoPen', tune.method="AIC",
ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5)
  
GLMvarISISscad(x, y, nsis=NULL, family=binomial(), folds=folds, 
rank.method="obj", tune.method="AIC", vartype="First", eps0=1e-5,
inittype='NoPen',  ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5)
}
\arguments{
  \item{x}{
     an (n * p)  matrix of features.
  }
  \item{y}{
     an (n) vector of response.
  }
  \item{nsis}{
  number of pedictors recuited by (I)SIS.
  }
  \item{family}{
    a description of the error distribution and link function to be used in the model.
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
   \code{ISIStypeCumulative} specifies whether to penalize variables selected by the previous step of the ISIS iteration in the following SCAD step. (\code{ISIStypeCumulative}= FALSE put penalties on all variables. In this case, the procedure is more likely to delete variables from previous step. )
   }
   \item{DOISIS}{
   \code{DOISIS} specifies whether to do iterative SIS.
   }
   \item{maxloop}{
   maximum number of loops in iterative SIS.
   }
  \item{eps0}{
      an effective zero relative to the scale of the maximum absolute marginal regression coefficients.
  }
}
\value{Returns an object with
\item{initRANKorder}{
    initial predictor ranking order for vanilla SIS.
    }
\item{detail.pickind, detail.ISISind}{
   details of each loop of ISIS.
   }
\item{normal.exit}{
   indicator of normal exit.
   }
\item{SISind}{
 the vector of indices selected by SIS.
  }
\item{ISISind}{
 the vector of indices selected by ISIS.
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

binom.van.sis=GLMvanISISscad(x, y, nsis, family=binomial(), tune.method='BIC')
binom.var.sis=GLMvarISISscad(x, y, nsis, family=binomial(), vartype='Second', 
tune.method='BIC')

#####compare the result
binom.van.sis$SIS

binom.van.sis$ISIS

binom.var.sis$SIS

binom.var.sis$ISIS
}

\keyword{file}
