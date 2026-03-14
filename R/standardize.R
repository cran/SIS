#' Standardization of High-Dimensional Design Matrices
#'
#' Standardizes the columns of a high-dimensional design matrix to mean zero
#' and unit Euclidean norm.
#'
#' Performs a location and scale transform to the columns of the original
#' design matrix, so that the resulting design matrix with \eqn{p}-dimensional
#' observations \eqn{\{x_i : i=1,...,n\}} of the form
#' \eqn{x_i=(x_{i1},x_{i2},...,x_{ip})} satisfies \eqn{\sum_{i=1}^{n} x_{ij} =
#' 0} and \eqn{\sum_{i=1}^{n} x_{ij}^{2} = 1} for \eqn{j=1,...,p}.
#'
#' @export
#' @param X A design matrix to be standardized.
#' @return A design matrix with standardized predictors or columns.
#' @author Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, and
#' Yichao Wu
#' @references
#' Diego Franco Saldana and Yang Feng (2018) SIS: An R package for Sure Independence Screening in
#' Ultrahigh Dimensional Statistical Models, \emph{Journal of Statistical Software}, \bold{83}, 2, 1-25.
#' @keywords models
#' @examples
#'\dontrun{
#' set.seed(0)
#' n <- 400
#' p <- 50
#' rho <- 0.5
#' corrmat <- diag(rep(1 - rho, p)) + matrix(rho, p, p)
#' corrmat[, 4] <- sqrt(rho)
#' corrmat[4, ] <- sqrt(rho)
#' corrmat[4, 4] <- 1
#' corrmat[, 5] <- 0
#' corrmat[5, ] <- 0
#' corrmat[5, 5] <- 1
#' cholmat <- chol(corrmat)
#' x <- matrix(rnorm(n * p, mean = 15, sd = 9), n, p)
#' x <- x %*% cholmat
#'
#' x.standard <- standardize(x)
#' }
standardize <- function(X) {
  center <- colMeans(X)
  X.c <- sweep(X, 2, center)
  unit.var <- sqrt(apply(X.c, 2, crossprod))
  val <- sweep(X.c, 2, unit.var, "/")
  return(val)
}
