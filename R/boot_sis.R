#' Calculates confidence intervals by using the quantile bootstrap approach. 
#' 
#' The algorithm first conducts 10 repetitions of a 200 iterations
#' bootstrap. Then, it checks whether the standard error of the mean of those 10 estimations for both upper and lower confidence intervals is
#' lower than the 5 % of the length of the interval for each variable. If this condition is met for more than 95 % of the variables,
#' the bootstrap sample is adequate and we use the constructed sample to calculate upper and lower confidence intervals. If the condition is not met for more than 5 % of the variables,
#' then the variability is too high and we need to increase the number of bootstrap repetitions. Thus, the process of testing 200 bootstrap samples 10 times is repeated and added to the previous bootstrap estimates.
#' This process is repeated until the condition is met for more than 95 % of the variables.
#'
#' @importFrom ncvreg ncvsurv
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom gcdnet cv.gcdnet
#' @importFrom gcdnet gcdnet
#' @importFrom glmnet cv.glmnet
#' @importFrom msaenet aenet
#' @import doParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom methods as
#' @importFrom stats sd
#' @export
#' @param x The design matrix, of dimensions n * p, without an intercept. Each row is an observation vector.
#' @param y The response vector of dimension n * 1. Quantitative for
#' \code{family='gaussian'}, non-negative counts for \code{family='poisson'},
#' binary (0-1) for \code{family='binomial'}. For \code{family='cox'}, \code{y}
#' should be an object of class \code{Surv}, as provided by the function
#' \code{Surv()} in the package \pkg{survival}.
#' @param family Response type (see above).
#' @param penalty The penalty to be applied in the regularized likelihood
#' subproblems. 'SCAD' (the default), 'MCP', 'lasso', 'enet' (elastic-net), 'msaenet' (multi-step adaptive elastic-net) and
#' 'aenet' (adaptive elastic-net) are provided.
#' @param sig significance threshold for the confidence interval
#' @param covars factor covariate names
#' @param probs Quantiles to compare for the effect estimation. By default quantiles are 10th and 90th
#' @param parallel Specifies whether to conduct parallel computing. If TRUE, it uses the parameter \code{parallel} from the \pkg{glmnet} package
#' to parallelize the computation
#' @return A data frame with columns \code{coef}, \code{CI_low},
#' \code{CI_up}, \code{Est}, \code{CI_low_perc}, and \code{CI_up_perc}.
#' @author Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, Arce Domingo-Relloso and Yichao Wu
#' @references Jerome Friedman and Trevor Hastie and Rob Tibshirani (2010)
#' Regularization Paths for Generalized Linear Models Via Coordinate Descent.
#' \emph{Journal of Statistical Software}, \bold{33}(1), 1-22.
#'
#' Noah Simon and Jerome Friedman and Trevor Hastie and Rob Tibshirani (2011)
#' Regularization Paths for Cox's Proportional Hazards Model Via Coordinate
#' Descent. \emph{Journal of Statistical Software}, \bold{39}(5), 1-13.
#'
#' Patrick Breheny and Jian Huang (2011) Coordiante Descent Algorithms for
#' Nonconvex Penalized Regression, with Applications to Biological Feature
#' Selection. \emph{The Annals of Applied Statistics}, \bold{5}, 232-253.
#'
#' Hirotogu Akaike (1973) Information Theory and an Extension of the Maximum
#' Likelihood Principle. In \emph{Proceedings of the 2nd International
#' Symposium on Information Theory}, BN Petrov and F Csaki (eds.), 267-281.
#'
#' Gideon Schwarz (1978) Estimating the Dimension of a Model. \emph{The Annals
#' of Statistics}, \bold{6}, 461-464.
#'
#' Jiahua Chen and Zehua Chen (2008) Extended Bayesian Information Criteria for
#' Model Selection with Large Model Spaces. \emph{Biometrika}, \bold{95},
#' 759-771.

boot_sis <- function(x, y, family, penalty, sig=0.05, covars, probs=c(0.1,0.9), parallel=TRUE) {
  if (penalty=='lasso'){
    fit <- cv.glmnet(x, y, family=family, parallel=parallel)
    lambda <- fit$lambda.1se
    coef = coef(fit, s = lambda) 
  } else if (penalty=='enet' | penalty=='msaenet' | penalty=='aenet'){
    fit <- cv.glmnet(x, y, alpha=0.05, family=family, parallel=parallel)
    lambda <- fit$lambda.1se
    coef <- coef(fit, s = lambda)
    if(any(row.names(coef)=='(Intercept)')){
      coef.aenet <- coef(fit, s = lambda)[-1,1]
    } else{
      coef.aenet <- coef(fit, s = lambda)[,1]
    }
    
    if (penalty=='msaenet'){
      fit_aenet <- aenet(x, y, family = family, init = "ridge", alphas = 0.05,
                         rule = "lambda.1se", parallel = parallel, verbose = TRUE, penalty.factor.init=(pmax(abs(coef.aenet), .Machine$double.eps))^(-1))
      coef <- data.frame(coef(fit_aenet))
      row.names(coef) <- colnames(x)
    } else if (penalty=='aenet'){
      if (family=='cox'){
        fit = Coxnet(x, y, penalty = "Enet", alpha = 0.05, nlambda = 50, nfolds = 10,
                     inzero = FALSE, adaptive = c(TRUE, FALSE), aini=list(wbeta=(pmax(abs(as.numeric(coef.aenet)), .Machine$double.eps))^(-1)))
        coef = as(as.matrix(fit[['Beta']]), "dgCMatrix")
        row.names(coef) <- colnames(x)
        lambda = fit[['lambda.max']]
      } else {
        if (family=='gaussian'){
          fit_gcdnet <- cv.gcdnet(x, y, nfolds=10, lambda2=0.95, standardize=TRUE, method='ls')
        } else {
          fit_gcdnet <- cv.gcdnet(x, y, nfolds=10, lambda2=0.95, standardize=TRUE)
        }
        lambda = fit_gcdnet$lambda.1se
        coef = coef(fit_gcdnet, s = "lambda.1se")
        if(any(row.names(coef)=='(Intercept)')){
          pf=as.vector(pmax(abs(coef), .Machine$double.eps)^(-1))[-1]
        } else{
          pf=as.vector(pmax(abs(coef), .Machine$double.eps)^(-1))
        }
      }
    }} else {
      if (family != 'cox') {
        fit = cv.ncvreg(x, y, family = family, penalty = penalty)
        cv.1se.ind = min(which(fit$cve<fit$cve[fit$min]+fit$cvse[fit$min]))
        coef = fit$fit$beta[, cv.1se.ind]  # extract coefficients at a single value of lambda, including the intercept
        coef <- data.frame(coef)
      } else {
        fit = ncvreg::cv.ncvsurv(x, y, family = family, penalty = penalty)
        cv.1se.ind = min(which(fit$cve<fit$cve[fit$min]+fit$cvse[fit$min]))
        coef = fit$fit$beta[, cv.1se.ind]
        coef <- data.frame(coef)
      }
    }
  
  
  ci_low <- ci_up <- bootstrap <- c()
  count <- 0
  k <- 0
  message('Getting bootstrap confidence intervals')
  
  repeat{
    for(j in 1:10){
      p <- ncol(x)
      n <- nrow(x)
      boot <- foreach (i=1:200, .combine='cbind') %dopar% {
        if(exists(".Random.seed")){
          rm(.Random.seed, envir=globalenv())
        }
        ind <- sample(1:n, replace=TRUE)
        if (penalty=='lasso'){
          fit.i <- glmnet(x[ind,], y[ind], family=family, lambda=lambda, parallel=parallel)
          b = coef(fit.i) 
          b
        } else if (penalty=='enet'){
          fit.i <- glmnet(x[ind,], y[ind], alpha=0.05, family=family, lambda=lambda, parallel=parallel)
          b = coef(fit.i)
          b
        } else if (penalty=='msaenet'){
          fit.i = aenet(x[ind,], y[ind], family = family, init = "ridge", alphas = 0.05,
                        rule = "lambda.1se", parallel = parallel, verbose = TRUE, penalty.factor.init=(pmax(abs(coef.aenet), .Machine$double.eps))^(-1))
          b = fit.i[['beta']]
          b
        } else if (penalty=='aenet'){
          if (family=='cox'){
            fit.i <- try(Coxnet(x[ind,], y[ind], penalty = "Enet", alpha = 0.05, lambda = lambda, nfolds = 10,
                                inzero = FALSE, adaptive = c(TRUE, FALSE), aini=list(wbeta=(pmax(abs(as.numeric(coef.aenet)), .Machine$double.eps))^(-1))))
            while(class(fit.i)=="try-error"){
              print('try')
              ind <- sample(1:n, replace=TRUE)
              fit.i <- try(Coxnet(x[ind,], y[ind], penalty = "Enet", alpha = 0.05, lambda = lambda, nfolds = 10,
                                  inzero = FALSE, adaptive = c(TRUE, FALSE), aini=list(wbeta=(pmax(abs(as.numeric(coef.aenet)), .Machine$double.eps))^(-1))))
            }
            b = as(as.matrix(fit.i[['Beta']]), "dgCMatrix")
            b
          } else {
            if (family=='gaussian'){
              fit.i <- gcdnet(x[ind,], y[ind], lambda=lambda, lambda2=0.95, pf=pf, standardize=TRUE, method='ls')
            } else {
              fit.i <- gcdnet(x[ind,], y[ind], lambda=lambda, lambda2=0.95, pf=pf, standardize=TRUE)
            }
            b = coef(fit.i)
            b
          }} else{
            if (family != 'cox') {
              fit = cv.ncvreg(x[ind,], y[ind], family = family, penalty = penalty)
              cv.1se.ind = min(which(fit$cve<fit$cve[fit$min]+fit$cvse[fit$min]))
              b = fit$fit$beta[, cv.1se.ind]  # extract coefficients at a single value of lambda, including the intercept
              b
            } else {
              fit = ncvreg::cv.ncvsurv(x[ind,], y[ind], family = family, penalty = penalty)
              cv.1se.ind = min(which(fit$cve<fit$cve[fit$min]+fit$cvse[fit$min]))
              b = fit$fit$beta[, cv.1se.ind]  # extract coefficients at a single value of lambda
              b
            }
          }
      }
      
      out <- as.data.frame(t(apply(boot, 1, quantile, probs=c(sig/2, 1-sig/2))))
      out[which(is.na(out[,1])),1] <- 0
      out[which(is.na(out[,2])),2] <- 0
      colnames(out) <- c('Lower', 'Upper')
      ci_low <- cbind(ci_low, out[,1])
      ci_up <- cbind(ci_up, out[,2])
      bootstrap <- cbind(bootstrap, boot)
      count <- count + 1
    }
    k <- k + 10
    message('Number of bootstrap repetition block:')
    print(k)
    if (length(which((apply(ci_low, 1, function(x) sd(x)/sqrt(k)) > 0.05*apply(abs(as.matrix(ci_up) - as.matrix(ci_low)), 1, mean)) | (apply(ci_up, 1, function(x) sd(x)/sqrt(k)) > 0.05*apply(abs(as.matrix(ci_low) - as.matrix(ci_up)), 1, mean))))<0.05*dim(x)[2]){
      break
    }}
  
  out <- as.data.frame(t(apply(bootstrap, 1, quantile, probs=c(sig/2, 1-sig/2))))
  if(!is.null(colnames(x))){
    if(any(rownames(out)=='(Intercept)')){
      out <- cbind(c('(Intercept)', colnames(x)), out)
    } else {
      out <- cbind(colnames(x), out)
    }
  } else{
    out <- cbind(row.names(out), out)
  }
  colnames(out) <- c('var', 'CI_low', 'CI_up')
  coef <- data.frame(rownames(coef), coef[,1])
  colnames(coef) <- c('var', 'coef')
  out <- merge(coef, out, by='var')
  
  if(family=='cox' | family=='binomial'){
    percentiles <- cbind(rep(NA, dim(out)[1]), rep(NA, dim(out)[1]))
    percentiles[which(as.character(out$var)!='(Intercept)'),] <- t(apply(x[,as.character(out$var)[which(as.character(out$var)!='(Intercept)')]], 2, function(x) quantile(x, probs=probs)))
    out$CI_up_perc <- out$CI_low_perc <- out$Est <- rep(NA, dim(out)[1])
    out$Est[which(as.character(out$var)!='(Intercept)')] <- exp((percentiles[which(as.character(out$var)!='(Intercept)'),2] - percentiles[which(as.character(out$var)!='(Intercept)'),1])*as.numeric(as.character(out$coef[!(as.character(out$var)=='(Intercept)')])))
    out$CI_low_perc[which(as.character(out$var)!='(Intercept)')] <- exp((percentiles[which(as.character(out$var)!='(Intercept)'),2] - percentiles[which(as.character(out$var)!='(Intercept)'),1])*as.numeric(as.character(out$CI_low[which(as.character(out$var)!='(Intercept)')])))
    out$CI_up_perc[which(as.character(out$var)!='(Intercept)')] <- exp((percentiles[which(as.character(out$var)!='(Intercept)'),2] - percentiles[which(as.character(out$var)!='(Intercept)'),1])*as.numeric(as.character(out$CI_up[which(as.character(out$var)!='(Intercept)')])))
    if(!is.null(covars)){
      for (i in covars){
        out[which(as.character(out$var)==i),'Est'] <- exp(as.numeric(out[which(as.character(out$var)==i), 'coef']))
        out[which(as.character(out$var)==i),'CI_low_perc'] <- exp(as.numeric(out[which(as.character(out$var)==i), 'CI_low']))
        out[which(as.character(out$var)==i),'CI_up_perc'] <- exp(as.numeric(out[which(as.character(out$var)==i), 'CI_up']))
      }
    }
  } else if (family=='gaussian'){
    percentiles <- cbind(rep(NA, dim(out)[1]), rep(NA, dim(out)[1]))
    percentiles[which(as.character(out$var)!='(Intercept)'),] <- t(apply(x[,as.character(out$var)[which(as.character(out$var)!='(Intercept)')]], 2, function(x) quantile(x, probs=probs)))
    out$CI_up_perc <- out$CI_low_perc <- out$Est <- rep(NA, dim(out)[1])
    out$Est[which(as.character(out$var)!='(Intercept)')] <- (percentiles[which(as.character(out$var)!='(Intercept)'),2] - percentiles[which(as.character(out$var)!='(Intercept)'),1])*as.numeric(as.character(out$coef[!(as.character(out$var)=='(Intercept)')]))
    out$CI_low_perc[which(as.character(out$var)!='(Intercept)')] <- (percentiles[which(as.character(out$var)!='(Intercept)'),2] - percentiles[which(as.character(out$var)!='(Intercept)'),1])*as.numeric(as.character(out$CI_low[which(as.character(out$var)!='(Intercept)')]))
    out$CI_up_perc[which(as.character(out$var)!='(Intercept)')] <- (percentiles[which(as.character(out$var)!='(Intercept)'),2] - percentiles[which(as.character(out$var)!='(Intercept)'),1])*as.numeric(as.character(out$CI_up[which(as.character(out$var)!='(Intercept)')]))
    if(!is.null(covars)){
      for (i in covars){
        out[which(as.character(out$var)==i),'Est'] <- as.numeric(out[which(as.character(out$var)==i), 'coef'])
        out[which(as.character(out$var)==i),'CI_low_perc'] <- as.numeric(out[which(as.character(out$var)==i), 'CI_low'])
        out[which(as.character(out$var)==i),'CI_up_perc'] <- as.numeric(out[which(as.character(out$var)==i), 'CI_up'])
      }
    }
  }
  return(cis=out)
}
