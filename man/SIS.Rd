\name{SIS}
\alias{SIS}
\title{
(Iterative) Sure Independence Screening ((I)SIS) and fitting in Generalized Linear Models and Cox 
Proportional Hazards Regression Models
}
\description{
 This function first implements the Iterative Sure Independence Screening for different variants of 
 (I)SIS, and then fits the final regression model using the R packages \pkg{ncvreg} and \pkg{glmnet}
 for the SCAD/MCP/LASSO regularized loglikelihood for the variables picked by (I)SIS.
}
\usage{
SIS(x, y, family = c("gaussian","binomial","poisson","cox"), 
penalty=c("SCAD","MCP","lasso"), concavity.parameter = 
switch(penalty, SCAD=3.7, 3), tune = c("cv","aic","bic","ebic"), 
nfolds = 10, type.measure = c("deviance","class","auc","mse",
"mae"), gamma.ebic = 1, nsis = NULL, iter = TRUE, iter.max = 
ifelse(greedy==FALSE,10,floor(nrow(x)/log(nrow(x)))), varISIS = 
c("vanilla","aggr","cons"), perm = FALSE, q = 1, greedy = FALSE, 
greedy.size = 1, seed = 0, standardize = TRUE)
}
\arguments{
  \item{x}{
     The design matrix, of dimensions n * p, without an intercept. Each row is an observation vector.
     \code{SIS} standardizes the data and includes an intercept by default.
  }
  \item{y}{
     The response vector of dimension n * 1. Quantitative for \code{family="gaussian"}, non-negative 
     counts for \code{family="poisson"}, binary (0-1) for \code{family=
     "binomial"}. For \code{family="cox"}, \code{y} should be an object of class \code{Surv}, 
     as provided by the function \code{Surv()} in the package \pkg{survival}.
  }
  \item{family}{
    Response type (see above).
  }
  \item{penalty}{
     The penalty to be applied in the regularized likelihood subproblems. "SCAD" (the default), "MCP", 
     or "lasso" are provided.
  }
  \item{concavity.parameter}{
     The tuning parameter used to adjust the concavity of the SCAD/MCP penalty. Default is 3.7 for 
     SCAD and 3 for MCP.
   }
  \item{tune}{
     Method for tuning the regularization parameter of the penalized likelihood subproblems and of the 
     final model selected by (I)SIS. Options include \code{tune="cv"}, \code{tune="aic"}, \code{tune="bic"}, 
     and \code{tune="ebic"}.
  }
  \item{nfolds}{
     Number of folds used in cross-validation. The default is 10.
  }
  \item{type.measure}{
      Loss to use for cross-validation. Currently five options, not all available for all models. The 
      default is \code{type.measure="deviance"}, which uses squared-error for gaussian models (also 
      equivalent to \code{type.measure="mse"} in this case), deviance for logistic and poisson regression, 
      and partial-likelihood for the Cox model. Both \code{type.measure="class"} and \code{type.measure="auc"} 
      apply only to logistic regression and give misclassification error and area under the ROC curve, 
      respectively. \code{type.measure="mse"} or \code{type.measure=
      "mae"} (mean absolute error) can be used by all models except the \code{"cox"}; they measure the 
      deviation from the fitted mean to the response. For \code{penalty=
      "SCAD"} and \code{penalty="MCP"}, only \code{type.measure="deviance"} is available.
  }
  \item{gamma.ebic}{
      Specifies the parameter in the Extended BIC criterion penalizing the size of the corresponding model 
      space. The default is \code{gamma.ebic=1}. See references at the end for details.
   }
  \item{nsis}{
      Number of pedictors recuited by (I)SIS.
  }
  \item{iter}{
    Specifies whether to perform iterative SIS. The default is \code{iter=TRUE}.
  }
  \item{iter.max}{
    Maximum number of iterations for (I)SIS and its variants.
  }
  \item{varISIS}{
    Specifies whether to perform any of the two ISIS variants based on randomly splitting the sample into 
    two groups. The variant \code{varISIS="aggr"} is an aggressive variable screening procedure, while 
    \code{varISIS="cons"} is a more conservative approach. The default is \code{varISIS="vanilla"}, which 
    performs the traditional vanilla version of ISIS. See references at the end for details.
  }
  \item{perm}{
   Specifies whether to impose a data-driven threshold in the size of the active sets calculated during 
   the ISIS procedures. The threshold is calculated by first decoupling the predictors \eqn{x_i} and 
   response \eqn{y_i} through a random permutation \eqn{\pi} of \eqn{(1,...,n)} to form a null model. 
   For this newly permuted data, marginal regression coefficients for each predictor are recalculated. 
   As the marginal regression coeffcients of the original data should be larger than most recalculated 
   coefficients in the null model, the data-driven threshold is given by the \eqn{q}th quantile of the 
   null coefficients. This data-driven threshold only allows a \eqn{1-q} proportion of inactive variables 
   to enter the model when \eqn{x_i} and \eqn{y_i} are not related (in the null model). The default is 
   here is \code{perm=FALSE}. See references at the end for details.
   }
  \item{q}{
   Quantile for calculating the data-driven threshold in the permutation-based ISIS. The default is 
   \code{q=1} (i.e., the maximum absolute value of the permuted estimates).
   }
   \item{greedy}{
   Specifies whether to run the greedy modification of the permutation-based ISIS. The default is 
   \code{greedy=FALSE}.
   }
   \item{greedy.size}{
   Maximum size of the active sets in the greedy modification of the permutation-based ISIS. The 
   default is \code{greedy.size=1}.
   }
   \item{seed}{
   Random seed used for sample splitting, random permutation, and cross-
   validation sampling of training and test sets.
   }
   \item{standardize}{
   Logical flag for x variable standardization, prior to performing (iterative) variable screening. 
   The resulting coefficients are always returned on the original scale. Default is \code{standardize=TRUE}. 
   If variables are in the same units already, you might not wish to standardize.
   }   
}
\value{Returns an object with
    \item{ix}{
    The vector of indices selected by (I)SIS.
    }   
    \item{coef.est}{
    The vector of coefficients of the final model selected by (I)SIS.
     } 
    \item{fit}{
    A fitted object of type \code{ncvreg}, \code{cv.ncvreg}, \code{glmnet}, or \code{cv.glmnet} for
    the final model selected by the (I)SIS procedure. If \code{tune="cv"}, the returned fitted object 
    is of type \code{cv.ncvreg} if \code{penalty="SCAD"} or \code{penalty="MCP"}; 
    otherwise, the returned fitted object is of type \code{cv.glmnet}. For the remaining 
    options of \code{tune}, the returned object is of type \code{glmnet} if \code{penalty="lasso"}, 
    and \code{ncvreg} otherwise.
    }
    \item{path.index}{
    The index along the solution path of \code{fit} for which the criterion specified in 
    \code{tune} is minimized.
    }
}
\references{
Jianqing Fan and Jinchi Lv (2008) Sure Independence Screening for Ultrahigh Dimensional Feature Space
(with discussion). \emph{Journal of Royal Statistical Society B}, \bold{70}, 849-911.

Jianqing Fan and Rui Song (2010) Sure Independence Screening in Generalized Linear Models with NP-Dimensionality. 
\emph{Journal of Royal Statistical Society B}, \bold{38}, 3567-3604.

Jianqing Fan, Richard Samworth, and Yichao Wu (2009) Ultrahigh Dimensional Feature Selection: Beyond the Linear
Model. \emph{Journal of Machine Learning Research}, \bold{10}, 2013-2038.

Jianqing Fan, Yang Feng, and Yichao Wu (2010) High-dimensional Variable Selection for Cox Proportional Hazards Model.
\emph{IMS Collections}, \bold{6}, 70-86.

Jianqing Fan, Yang Feng, and Rui Song (2011) Nonparametric Independence Screening in Sparse Ultrahigh Dimensional Additive Models.
\emph{Journal of the American Statistical Association}, \bold{106}, 544-557.

Diego Franco Saldana and Yang Feng (2014) Sure Independence Screening in Ultrahigh Dimensional Statistical Models, manuscript.

Jiahua Chen and Zehua Chen (2008) Extended Bayesian Information Criteria for Model Selection with Large Model Spaces.
\emph{Biometrika}, \bold{95}, 759-771.
}
\author{Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, and Yichao Wu}
\seealso{
\code{\link{predict.SIS}}, \code{\link{print.SIS}}
}
\examples{

set.seed(0)
n = 400; p = 50; rho = 0.5
corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
corrmat[,4] = sqrt(rho)
corrmat[4, ] = sqrt(rho)
corrmat[4,4] = 1
corrmat[,5] = 0
corrmat[5, ] = 0
corrmat[5,5] = 1
cholmat = chol(corrmat)
x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
x = x\%*\%cholmat

# gaussian response 
set.seed(1)
b = c(4,4,4,-6*sqrt(2),4/3)
y=x[, 1:5]\%*\%b + rnorm(n)
model11=SIS(x, y, family="gaussian", tune="bic")
model12=SIS(x, y, family="gaussian", tune="bic", varISIS="aggr", seed=11)
model11$ix
model12$ix

# binary response 
set.seed(2)
feta = x[, 1:5]\%*\%b; fprob = exp(feta)/(1+exp(feta))
y = rbinom(n, 1, fprob)
model21=SIS(x, y, family="binomial", tune="bic")
model22=SIS(x, y, family="binomial", tune="bic", varISIS="aggr", seed=21)
model21$ix
model22$ix

# poisson response
set.seed(3)
b = c(0.6,0.6,0.6,-0.9*sqrt(2))
myrates = exp(x[, 1:4]\%*\%b)
y = rpois(n, myrates)
model31=SIS(x, y, family="poisson", tune="bic", perm=TRUE, q=0.9, 
            greedy=TRUE, seed=31)
model32=SIS(x, y, family="poisson", tune="bic", varISIS="aggr", 
            perm=TRUE, q=0.9, seed=32)
model31$ix
model32$ix

# Cox model
set.seed(4)
b = c(4,4,4,-6*sqrt(2),4/3)
myrates = exp(x[, 1:5]\%*\%b)
Sur = rexp(n,myrates); CT = rexp(n,0.1)
Z = pmin(Sur,CT); ind = as.numeric(Sur<=CT)
y = Surv(Z,ind)
model41=SIS(x, y, family="cox", penalty="lasso", tune="bic", 
            varISIS="aggr", seed=41)
model42=SIS(x, y, family="cox", penalty="lasso", tune="bic", 
            varISIS="cons", seed=41)
model41$ix
model42$ix

}
\keyword{file}

