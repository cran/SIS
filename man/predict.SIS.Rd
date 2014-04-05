\name{predict.SIS}
\alias{predict.SIS}
\title{
  Model prediction based on a fitted SIS object.
}
\description{
Similar to the usual predict methods, this function returns predictions from a fitted \code{"SIS"} object.
}
\usage{
\method{predict}{SIS}(object, newx, s=object$path.index, type = c("response", 
"link", "class"), ...)
}
\arguments{
  \item{object}{
     Fitted \code{"SIS"} model object.
}
  \item{newx}{
     Matrix of new values for \code{x} at which predictions are to be made, without the intercept term.
}
  \item{s}{
     Indices of the penalty parameter \code{lambda} of the final fitted model by (I)SIS at which predictions are required. By default, 
     only the index minimizing the criterion \code{tune} is returned.
}
  \item{type}{
    Type of prediction required. Type \code{"response"} gives the fitted values for \code{"gaussian"}, fitted probabilities for
    \code{"binomial"}, fitted mean for 
    \code{"poisson"}, and the fitted relative risk for \code{"cox"}. Type \code{"link"} returns the linear predictors for 
    \code{"binomial"}, \code{"poisson"} and \code{"cox"} models; for \code{"gaussian"} models it is equivalent to type 
    \code{"response"}. Type \code{"class"} applies only to \code{"binomial"} models, and produces the class label corresponding 
    to the maximum probability (0-1 labels).
    }
  \item{\dots}{Not used. Other arguments to predict.}
}
\value{  
The object returned depends on type.
}
\author{Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, and Yichao Wu}
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
}
\seealso{
\code{\link{SIS}}, \code{\link{print.SIS}}
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
testX = matrix(rnorm(10*p, mean=0, sd=1), nrow=10, ncol=p)

# gaussian response 
set.seed(1)
b = c(4,4,4,-6*sqrt(2),4/3)
y=x[, 1:5]\%*\%b + rnorm(n)
model1=SIS(x, y, family="gaussian", tune="bic", varISIS="aggr", seed=11)

predict(model1, testX, type="response")
predict(model1, testX, s=1:10, type="response")

# binary response 
set.seed(2)
feta = x[, 1:5]\%*\%b; fprob = exp(feta)/(1+exp(feta))
y = rbinom(n, 1, fprob)
model2=SIS(x, y, family="binomial", tune="bic", varISIS="aggr", seed=21)

predict(model2, testX, type="response")
predict(model2, testX, type="link")
predict(model2, testX, type="class")

predict(model2, testX, s=1:10, type="response")
predict(model2, testX, s=1:10, type="link")
predict(model2, testX, s=1:10, type="class")

# poisson response
set.seed(3)
b = c(0.6,0.6,0.6,-0.9*sqrt(2))
myrates = exp(x[, 1:4]\%*\%b)
y = rpois(n, myrates)
model3=SIS(x, y, family="poisson", tune="bic", varISIS="aggr", seed=31)

predict(model3, testX, type="response")
predict(model3, testX, type="link")

predict(model3, testX, s=1:10, type="response")
predict(model3, testX, s=1:10, type="link")

# Cox model
set.seed(4)
b = c(4,4,4,-6*sqrt(2),4/3)
myrates = exp(x[, 1:5]\%*\%b)
Sur = rexp(n,myrates); CT = rexp(n,0.1)
Z = pmin(Sur,CT); ind = as.numeric(Sur<=CT)
y = Surv(Z,ind)
model4=SIS(x, y, family="cox", penalty="lasso", tune="bic", 
           varISIS="aggr", seed=41)

predict(model4, testX, type="response")
predict(model4, testX, type="link")

predict(model4, testX, s=1:10, type="response")
predict(model4, testX, s=1:10, type="link")

}
\keyword{file}
