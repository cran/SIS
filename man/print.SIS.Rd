\name{print.SIS}
\alias{print.SIS}
\title{Print a SIS object}
\description{
Print a summary of the selected predictors and its estimated coefficients.
}
\usage{
\method{print}{SIS}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{Fitted \code{"SIS"} model object.}
  \item{digits}{Significant digits in printout.}
  \item{\dots}{Additional print arguments.}
}
\details{
  The call that produced object \code{x} is printed, followed by the indices of the selected predictors (i.e., \code{x$ix})
  and its estimated coefficients (i.e., \code{x$coef.est}) for the final model picked up by (I)SIS.
}
\value{
  The fitted regression object of type \code{"cv.ncvreg"}, \code{"ncvreg"}, \code{"cv.glmnet"}, or \code{"glmnet"} and
  its corresponding path index where the minimum criterion in \code{tune} occurs are silently returned.
}
\author{Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, and Yichao Wu}
\references{
Jianqing Fan and Jinchi Lv (2008) Sure Independence Screening for Ultrahigh Dimensional Feature Space
(with discussion). \emph{Journal of Royal Statistical Society B}, \bold{70}, 849-911.

Jianqing Fan and Rui Song (2010) Sure Independence Screening in Generalized Linear Models with NP-Dimensionality. 
\emph{The Annals of Statistics}, \bold{38}, 3567-3604.

Jianqing Fan, Richard Samworth, and Yichao Wu (2009) Ultrahigh Dimensional Feature Selection: Beyond the Linear
Model. \emph{Journal of Machine Learning Research}, \bold{10}, 2013-2038.

Jianqing Fan, Yang Feng, and Yichao Wu (2010) High-dimensional Variable Selection for Cox Proportional Hazards Model.
\emph{IMS Collections}, \bold{6}, 70-86.

Jianqing Fan, Yang Feng, and Rui Song (2011) Nonparametric Independence Screening in Sparse Ultrahigh Dimensional Additive Models.
\emph{Journal of the American Statistical Association}, \bold{106}, 544-557.

Diego Franco Saldana and Yang Feng (2014) Sure Independence Screening in Ultrahigh Dimensional Statistical Models, manuscript.
}
\seealso{
\code{\link{SIS}}, \code{\link{predict.SIS}}
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
model1=SIS(x, y, family="gaussian", tune="bic", varISIS="aggr", seed=11)

print(model1, digits=3)
print(model1, digits=7)
}
\keyword{models}
