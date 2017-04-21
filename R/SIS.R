#' (Iterative) Sure Independence Screening ((I)SIS) and Fitting in Generalized
#' Linear Models and Cox's Proportional Hazards Models
#' 
#' This function first implements the Iterative Sure Independence Screening for
#' different variants of (I)SIS, and then fits the final regression model using
#' the R packages \pkg{ncvreg} and \pkg{glmnet} for the SCAD/MCP/LASSO
#' regularized loglikelihood for the variables picked by (I)SIS.
#' 
#' @export
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#' @importFrom glmnet cv.glmnet
#' @importFrom ncvreg cv.ncvreg
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom stats cor glm.fit
#' @importFrom stats glm.fit
#' @importFrom stats quantile
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats gaussian
#' @importFrom stats poisson
#' @importFrom stats deviance
#' @importFrom stats predict
#' 
#' 
#' @param x The design matrix, of dimensions n * p, without an intercept. Each
#' row is an observation vector.  \code{SIS} standardizes the data and includes
#' an intercept by default.
#' @param y The response vector of dimension n * 1. Quantitative for
#' \code{family='gaussian'}, non-negative counts for \code{family='poisson'},
#' binary (0-1) for \code{family='binomial'}. For \code{family='cox'}, \code{y}
#' should be an object of class \code{Surv}, as provided by the function
#' \code{Surv()} in the package \pkg{survival}.
#' @param family Response type (see above).
#' @param penalty The penalty to be applied in the regularized likelihood
#' subproblems. 'SCAD' (the default), 'MCP', or 'lasso' are provided.
#' @param concavity.parameter The tuning parameter used to adjust the concavity
#' of the SCAD/MCP penalty. Default is 3.7 for SCAD and 3 for MCP.
#' @param tune Method for tuning the regularization parameter of the penalized
#' likelihood subproblems and of the final model selected by (I)SIS. Options
#' include \code{tune='bic'}, \code{tune='ebic'}, \code{tune='aic'}, and
#' \code{tune='cv'}.
#' @param nfolds Number of folds used in cross-validation. The default is 10.
#' @param type.measure Loss to use for cross-validation. Currently five
#' options, not all available for all models. The default is
#' \code{type.measure='deviance'}, which uses squared-error for gaussian models
#' (also equivalent to \code{type.measure='mse'} in this case), deviance for
#' logistic and poisson regression, and partial-likelihood for the Cox model.
#' Both \code{type.measure='class'} and \code{type.measure='auc'} apply only to
#' logistic regression and give misclassification error and area under the ROC
#' curve, respectively. \code{type.measure='mse'} or \code{type.measure='mae'}
#' (mean absolute error) can be used by all models except the \code{'cox'};
#' they measure the deviation from the fitted mean to the response. For
#' \code{penalty='SCAD'} and \code{penalty='MCP'}, only
#' \code{type.measure='deviance'} is available.
#' @param gamma.ebic Specifies the parameter in the Extended BIC criterion
#' penalizing the size of the corresponding model space. The default is
#' \code{gamma.ebic=1}. See references at the end for details.
#' @param nsis Number of pedictors recuited by (I)SIS.
#' @param iter Specifies whether to perform iterative SIS. The default is
#' \code{iter=TRUE}.
#' @param iter.max Maximum number of iterations for (I)SIS and its variants.
#' @param varISIS Specifies whether to perform any of the two ISIS variants
#' based on randomly splitting the sample into two groups. The variant
#' \code{varISIS='aggr'} is an aggressive variable screening procedure, while
#' \code{varISIS='cons'} is a more conservative approach. The default is
#' \code{varISIS='vanilla'}, which performs the traditional vanilla version of
#' ISIS. See references at the end for details.
#' @param perm Specifies whether to impose a data-driven threshold in the size
#' of the active sets calculated during the ISIS procedures. The threshold is
#' calculated by first decoupling the predictors \eqn{x_i} and response
#' \eqn{y_i} through a random permutation \eqn{\pi} of \eqn{(1,...,n)} to form
#' a null model.  For this newly permuted data, marginal regression
#' coefficients for each predictor are recalculated.  As the marginal
#' regression coeffcients of the original data should be larger than most
#' recalculated coefficients in the null model, the data-driven threshold is
#' given by the \eqn{q}th quantile of the null coefficients. This data-driven
#' threshold only allows a \eqn{1-q} proportion of inactive variables to enter
#' the model when \eqn{x_i} and \eqn{y_i} are not related (in the null model).
#' The default is here is \code{perm=FALSE}. See references at the end for
#' details.
#' @param q Quantile for calculating the data-driven threshold in the
#' permutation-based ISIS. The default is \code{q=1} (i.e., the maximum
#' absolute value of the permuted estimates).
#' @param greedy Specifies whether to run the greedy modification of the
#' permutation-based ISIS. The default is \code{greedy=FALSE}.
#' @param greedy.size Maximum size of the active sets in the greedy
#' modification of the permutation-based ISIS. The default is
#' \code{greedy.size=1}.
#' @param seed Random seed used for sample splitting, random permutation, and
#' cross- validation sampling of training and test sets.
#' @param standardize Logical flag for x variable standardization, prior to
#' performing (iterative) variable screening.  The resulting coefficients are
#' always returned on the original scale. Default is \code{standardize=TRUE}.
#' If variables are in the same units already, you might not wish to
#' standardize.
#' @return Returns an object with \item{ix}{ The vector of indices selected by
#' (I)SIS.  } \item{coef.est}{ The vector of coefficients of the final model
#' selected by (I)SIS.  } \item{fit}{ A fitted object of type \code{ncvreg},
#' \code{cv.ncvreg}, \code{glmnet}, or \code{cv.glmnet} for the final model
#' selected by the (I)SIS procedure. If \code{tune='cv'}, the returned fitted
#' object is of type \code{cv.ncvreg} if \code{penalty='SCAD'} or
#' \code{penalty='MCP'}; otherwise, the returned fitted object is of type
#' \code{cv.glmnet}. For the remaining options of \code{tune}, the returned
#' object is of type \code{glmnet} if \code{penalty='lasso'}, and \code{ncvreg}
#' otherwise.  } \item{path.index}{ The index along the solution path of
#' \code{fit} for which the criterion specified in \code{tune} is minimized.  }
#' @author Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, and
#' Yichao Wu
#' @seealso \code{\link{predict.SIS}}
#' @references 
#' Diego Franco Saldana and Yang Feng (2016) SIS: An R package for Sure Independence Screening in
#' Ultrahigh Dimensional Statistical Models, \emph{Journal of Statistical Software}, to appear.
#' 
#' Jianqing Fan and Jinchi Lv (2008) Sure Independence Screening
#' for Ultrahigh Dimensional Feature Space (with discussion). \emph{Journal of
#' Royal Statistical Society B}, \bold{70}, 849-911.
#' 
#' Jianqing Fan and Rui Song (2010) Sure Independence Screening in Generalized
#' Linear Models with NP-Dimensionality.  \emph{The Annals of Statistics},
#' \bold{38}, 3567-3604.
#' 
#' Jianqing Fan, Richard Samworth, and Yichao Wu (2009) Ultrahigh Dimensional
#' Feature Selection: Beyond the Linear Model. \emph{Journal of Machine
#' Learning Research}, \bold{10}, 2013-2038.
#' 
#' Jianqing Fan, Yang Feng, and Yichao Wu (2010) High-dimensional Variable
#' Selection for Cox Proportional Hazards Model. \emph{IMS Collections},
#' \bold{6}, 70-86.
#' 
#' Jianqing Fan, Yang Feng, and Rui Song (2011) Nonparametric Independence
#' Screening in Sparse Ultrahigh Dimensional Additive Models. \emph{Journal of
#' the American Statistical Association}, \bold{106}, 544-557.
#' 
#' Jiahua Chen and Zehua Chen (2008) Extended Bayesian Information Criteria for
#' Model Selection with Large Model Spaces. \emph{Biometrika}, \bold{95},
#' 759-771.
#' @keywords models
#' @examples
#' 
#' 
#' set.seed(0)
#' n = 400; p = 50; rho = 0.5
#' corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
#' corrmat[,4] = sqrt(rho)
#' corrmat[4, ] = sqrt(rho)
#' corrmat[4,4] = 1
#' corrmat[,5] = 0
#' corrmat[5, ] = 0
#' corrmat[5,5] = 1
#' cholmat = chol(corrmat)
#' x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
#' x = x%*%cholmat
#' 
#' # gaussian response 
#' set.seed(1)
#' b = c(4,4,4,-6*sqrt(2),4/3)
#' y=x[, 1:5]%*%b + rnorm(n)
#' model11=SIS(x, y, family='gaussian', tune='bic')
#' model12=SIS(x, y, family='gaussian', tune='bic', varISIS='aggr', seed=11)
#' model11$ix
#' model12$ix
#' 
#' # binary response 
#' set.seed(2)
#' feta = x[, 1:5]%*%b; fprob = exp(feta)/(1+exp(feta))
#' y = rbinom(n, 1, fprob)
#' model21=SIS(x, y, family='binomial', tune='bic')
#' model22=SIS(x, y, family='binomial', tune='bic', varISIS='aggr', seed=21)
#' model21$ix
#' model22$ix
#' 
#' # poisson response
#' set.seed(3)
#' b = c(0.6,0.6,0.6,-0.9*sqrt(2))
#' myrates = exp(x[, 1:4]%*%b)
#' y = rpois(n, myrates)
#' model31=SIS(x, y, family='poisson', penalty = 'lasso', tune='bic', perm=TRUE, q=0.9, 
#'             greedy=TRUE, seed=31)
#' #model32=SIS(x, y, family='poisson', penalty = 'lasso',  tune='bic', varISIS='aggr', 
#' #            perm=TRUE, q=0.9, seed=32)
#' model31$ix
#' #model32$ix
#' 
#' # Cox model
#' #set.seed(4)
#' #b = c(4,4,4,-6*sqrt(2),4/3)
#' #myrates = exp(x[, 1:5]%*%b)
#' #Sur = rexp(n,myrates); CT = rexp(n,0.1)
#' #Z = pmin(Sur,CT); ind = as.numeric(Sur<=CT)
#' #y = survival::Surv(Z,ind)
#' #model41=SIS(x, y, family='cox', penalty='lasso', tune='bic', 
#' #           varISIS='aggr', seed=41)
#' #model42=SIS(x, y, family='cox', penalty='lasso', tune='bic', 
#' #             varISIS='cons', seed=41)
#' #model41$ix
#' #model42$ix
#' 
#' 
SIS <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox"), penalty = c("SCAD", "MCP", "lasso"), 
    concavity.parameter = switch(penalty, SCAD = 3.7, 3), tune = c("bic", "ebic", "aic", "cv"), nfolds = 10, 
    type.measure = c("deviance", "class", "auc", "mse", "mae"), gamma.ebic = 1, nsis = NULL, iter = TRUE, iter.max = ifelse(greedy == 
        FALSE, 10, floor(nrow(x)/log(nrow(x)))), varISIS = c("vanilla", "aggr", "cons"), perm = FALSE, q = 1, 
    greedy = FALSE, greedy.size = 1, seed = 0, standardize = TRUE) {
    
    this.call = match.call()
    family = match.arg(family)
    penalty = match.arg(penalty)
    tune = match.arg(tune)
    type.measure = match.arg(type.measure)
    varISIS = match.arg(varISIS)
    
    if (is.null(x) || is.null(y)) 
        stop("The data is missing!")
    if (class(concavity.parameter) != "numeric") 
        stop("concavity.parameter must be numeric!")
    if (class(nfolds) != "numeric") 
        stop("nfolds must be numeric!")
    if (class(seed) != "numeric") 
        stop("seed must be numeric!")
    
    if (family == "cox" && penalty %in% c("SCAD", "MCP")) 
        stop("Cox model currently not implemented with selected penalty")
    
    if (type.measure %in% c("class", "auc") && family %in% c("gaussian", "poisson", "cox")) 
        stop("'class' and 'auc' type measures are only available for logistic regression")
    
    if (type.measure %in% c("class", "auc", "mse", "mae") && penalty %in% c("SCAD", "MCP")) 
        stop("Only 'deviance' is available as type.measure for non-convex penalties")
    
    fit = switch(family, gaussian = sisglm(x, y, "gaussian", penalty, concavity.parameter, tune, nfolds, type.measure, 
        gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize), binomial = sisglm(x, 
        y, "binomial", penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, 
        varISIS, perm, q, greedy, greedy.size, seed, standardize), poisson = sisglm(x, y, "poisson", penalty, 
        concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, 
        greedy, greedy.size, seed, standardize), cox = sisglm(x, y, "cox", penalty, concavity.parameter, tune, 
        nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, 
        standardize))
    fit$call = this.call
    class(fit) = c(class(fit), "SIS")
    return(fit)
}

sisglm <- function(x, y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, 
    iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1 = NULL, s2 = NULL, split.tries = 0) {
    
    storage.mode(x) = "numeric"
    n = dim(x)[1]
    p = dim(x)[2]
    models = vector("list")
    if (is.null(nsis) == TRUE) 
        nsis = calculate.nsis(family = family, varISIS = varISIS, n = n, p = p)
    if (is.null(s1) == TRUE) {
        set.seed(seed)
        split.sample = sample(1:n)
        s1 = split.sample[1:ceiling(n/2)]
        s2 = setdiff(split.sample, s1)
    }
    old.x = x
    if (standardize == TRUE) {
        x = scale(x)
    }
    iterind = 0
    
    if (iter == TRUE) {
        ix0 = sort(obtain.ix0(x = x, y = y, s1 = s1, s2 = s2, family = family, nsis = nsis, iter = iter, varISIS = varISIS, 
            perm = perm, q = q, greedy = greedy, greedy.size = greedy.size, iterind = iterind))
        repeat {
            iterind = iterind + 1
            cat("Iter", iterind, ", screening: ", ix0, "\n")
            if(length(ix0) == 1 & penalty == 'lasso'){
              ix0 = c(ix0, p+1-ix0)
            }
            pen.ind = ix0
            selection.fit = tune.fit(old.x[,ix0,drop = FALSE], y, family , penalty , concavity.parameter, tune, nfolds , type.measure , gamma.ebic)
            coef.beta = selection.fit$beta
            a0 = selection.fit$a0
            
            lambda  = selection.fit$lambda
            lambda.ind = selection.fit$lambda.ind
            ix1 = sort(ix0[selection.fit$ix])
            if (length(ix1) == 0) {
                split.tries = split.tries + 1
                split.sample = sample(1:n)
                s1 = split.sample[1:ceiling(n/2)]
                s2 = setdiff(split.sample, s1)
                cat("Sample splitting attempt: ", split.tries, "\n")
                if (split.tries >= 20) {
                  cat("No variables remaining after ", split.tries, " sample splitting attempts! \n")
                  cat("You can try a more conservative variable screening approach! \n")
                } else return(sisglm(old.x, y, family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic,
                  nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1, s2, split.tries))
            }
          
            
            cat("Iter", iterind, ", selection: ", ix1, "\n")
            if (length(ix1) >= nsis || iterind >= iter.max) {
                ix0 = ix1
                if (length(ix1) >= nsis) 
                  cat("Maximum number of variables selected \n")
                if (iterind >= iter.max) 
                  cat("Maximum number of iterations reached \n")
                break
            }
            
            models[[iterind]] = ix1
            flag.models = 0
            if (iterind > 1) {
                for (j in 1:(iterind - 1)) {
                  if (identical(models[[j]], ix1) == TRUE) 
                    flag.models = 1
                }
            }
            if (flag.models == 1) {
                cat("Model already selected \n")
                break
            }
            
            candind = setdiff(1:p, ix1)
            pleft = nsis - length(ix1)
            newix = sort(obtain.newix(x = x, y = y, candind = candind, ix1 = ix1, s1 = s1, s2 = s2, family = family, 
                pleft = pleft, varISIS = varISIS, perm = perm, q = q, greedy = greedy, greedy.size = greedy.size, 
                iterind = iterind))
            cat("Iter", iterind, ", conditional-screening: ", newix, "\n")
            ix1 = sort(c(ix1, newix))
            if (setequal(ix1, ix0)) {
               flag.models = 1
            }
            ix0 = ix1
            if(length(ix1) == 0) break
        }  # end repeat
        
    } else {
        # end if(iter==TRUE)
        ix0 = sort(obtain.ix0(x = x, y = y, s1 = s1, s2 = s2, family = family, nsis = nsis, iter = iter, varISIS = varISIS, 
            perm = perm, q = q, greedy = greedy, greedy.size = greedy.size, iterind = iterind))
        if(length(ix0) == 1 & penalty == 'lasso'){
          ix0 = c(ix0, p+1-ix0)
        }
        pen.ind = ix0
        selection.fit = tune.fit(old.x[, ix0, drop = FALSE], y, family , penalty , concavity.parameter, tune, nfolds , type.measure , gamma.ebic)
        coef.beta = selection.fit$beta
        a0 = selection.fit$a0
        lambda = selection.fit$lambda
        lambda.ind = selection.fit$lambda.ind
        ix1 = sort(ix0[selection.fit$ix])
    }
    
    
    
    if (family == "cox") {
      if (length(ix1) > 0){
      names(coef.beta) = paste("X", ix1, sep = "") 
      }
    }  else {
      coef.beta = c(a0, coef.beta) 
      if(length(ix1)>0){
        names(coef.beta) = c("(Intercept)", paste("X", ix1, sep = ""))
      }
    }
    
    return(list(ix = ix1, coef.est = coef.beta, fit = selection.fit$fit, lambda = lambda, lambda.ind = lambda.ind, ix0 = pen.ind))
}
