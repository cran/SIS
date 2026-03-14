################################
#####  Cox model           #####
#####  Enet & Lasso        #####
#####  Coordinate Descent  #####
################################


##########################################################
#####  Cox with Regularization (coordinate descent)  #####
##########################################################
###  lambda1*||_1+lambda2/2*||_2 <-> lambda*(alpha_1+(1-alpha)/2_2)

#' @importFrom stats weighted.mean

Coxnet <- function(x, y, penalty = c("Lasso", "Enet"), alpha = 1, lambda = NULL,
                   nlambda = 50, rlambda = NULL, nfolds = 1, foldid = NULL,
                   inzero = TRUE, adaptive = c(FALSE, TRUE), aini = NULL,
                   isd = FALSE, ifast = TRUE, keep.beta = FALSE, thresh = 1e-6,
                   maxit = 1e+5) {
  penalty <- match.arg(penalty)
  if (penalty == "Lasso") {
    penalty <- "Enet"
    alpha <- 1
  }

  fit <- coxEnet(
    x, y, alpha, lambda, nlambda, rlambda, nfolds, foldid, inzero,
    adaptive[1], aini, isd, ifast, keep.beta, thresh, maxit
  )

  class(fit) <- "Coxnet"
  fit
}


##########################
#####  Enet (Lasso)  #####
##########################

coxEnet <- function(x, y, alpha = 1, lambda = NULL, nlambda = 100,
                    rlambda = NULL, nfolds = 1, foldid = NULL, itrunc = TRUE,
                    adaptive = FALSE, aini = NULL, isd = FALSE, ifast = TRUE,
                    keep.beta = FALSE, thresh = 1e-7, maxit = 1e+5) {
  penalty <- ifelse(alpha == 1, "Lasso", "Enet")

  N0 <- nrow(x)
  p <- ncol(x)
  ifast <- as.integer(ifast)

  tem <- scaleC(x)
  xscale <- tem$sd
  x <- tem$x
  rm(tem)

  prep0 <- coxprep(x, y)

  if (adaptive) {
    if (is.null(aini)) {
      aini <- coxini(x, y)
    }
    wbeta <- aini$wbeta
    rm(aini)
  } else {
    wbeta <- rep(1, p)
  }

  if (is.null(lambda)) {
    lambda_max <- max_lambdaC(
      prep0$x, prep0$tevent, prep0$N, prep0$nevent, prep0$nevent1,
      prep0$loc1, prep0$n, alpha, wbeta, N0
    )
    lambda_min <- ifelse(
      is.null(rlambda),
      ifelse(N0 > p, lambda_max * 0.0001, lambda_max * 0.01),
      lambda_max * rlambda
    )
    lambda <- lambda_max * (lambda_min / lambda_max) ^
      (c(0:(nlambda - 1)) / (nlambda - 1))
  } else {
    nlambda <- length(lambda)
  }

  out <- coxenetC(
    prep0$x, prep0$tevent, alpha, lambda, nlambda, wbeta, prep0$N,
    prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, p, N0, thresh,
    maxit, ifast
  )
  nlambdai <- out$nlambda
  if (nlambdai == 0) {
    return(NULL)
  }
  lambdai <- lambda[1:nlambdai]
  if (!isd) {
    out$Beta <- matrix(out$Beta[, 1:nlambdai] / xscale, ncol = nlambdai)
  }
  out$nzero <- apply(out$Beta != 0, 2, sum)
  out$flag <- out$flag[1:nlambdai]

  if (nfolds == 1 & is.null(foldid)) {
    fit <- data.frame(lambda = lambdai, nzero = out$nzero)
    return(list(
      Beta = out$Beta, fit = fit, penalty = penalty, adaptive = adaptive,
      flag = out$flag
    ))
  }

  if (is.null(foldid)) {
    foldid <- coxsplit(y, nfolds)
  } else {
    nfolds <- max(foldid)
  }
  tb <- table(foldid)
  N0i <- numeric(nfolds)
  for (i in 1:nfolds) {
    N0i[i] <- sum(tb[-i])
  }

  prepk <- list()
  for (i in 1:nfolds) {
    temid <- which(foldid != i)
    prepk[[i]] <- coxprep(x[temid, ], y[temid, ])
  }
  weighti <- as.vector(tapply(y[, "status"], foldid, sum))

  outi <- list()
  cvPL <- matrix(NA, nrow = nfolds, ncol = nlambdai)
  for (i in 1:nfolds) {
    outi[[i]] <- cvcoxenetC(
      prepk[[i]]$x, prepk[[i]]$tevent, alpha, lambdai, nlambdai, wbeta,
      prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1,
      prepk[[i]]$n, p, N0i[i], thresh, maxit, 0, prep0$x, prep0$N,
      prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n
    )
    cvPL[i, 1:outi[[i]]$nlambda] <- outi[[i]]$lf[1:outi[[i]]$nlambda] -
      outi[[i]]$ll[1:outi[[i]]$nlambda]
  }

  cvPL <- matrix(cvPL[, 1:nlambdai], ncol = nlambdai)
  cvraw <- cvPL / weighti
  nfoldi <- apply(!is.na(cvraw), 2, sum)
  rm(cvPL)
  cvm <- apply(cvraw, 2, weighted.mean, w = weighti, na.rm = TRUE)
  cvse <- sqrt(
    apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w = weighti,
          na.rm = TRUE) / (nfoldi - 1)
  )

  indexi <- which.max(cvm)
  temi <- rep("", nlambdai)
  temi[indexi] <- "max"
  temCV <- data.frame(
    lambda = lambdai, cvm = cvm, cvse = cvse, nzero = out$nzero,
    index = temi, stringsAsFactors = FALSE
  )

  if (!itrunc) {
    rm(outi)
    if (!keep.beta) {
      return(list(
        Beta = out$Beta[, indexi], fit = temCV, lambda.max = lambdai[indexi],
        penalty = penalty, adaptive = adaptive, flag = out$flag
      ))
    }
    return(list(
      Beta = out$Beta, fit = temCV, lambda.max = lambdai[indexi],
      penalty = penalty, adaptive = adaptive, flag = out$flag
    ))
  }

  il0 <- indexi
  cvm <- list()
  cv.max <- rep(NA, nlambdai)
  repeat {
    numi <- out$nzero[il0]
    Betai <- sapply(outi, function(x) x$Beta[, il0])
    Betao <- apply(Betai != 0, 2, sum)
    numi2 <- min(max(Betao), numi)

    if (numi2 > 0) {
      cvPL <- matrix(NA, nrow = nfolds, ncol = numi2)
      for (i in 1:nfolds) {
        Betaj <- Betai[, i]
        numj <- min(Betao[i], numi)
        if (numj == 0) {
          cvPL[i, ] <- cvtrimC(
            c(0.0, 0.0), numj, numi2, c(0, 0), prep0$x, prep0$N,
            prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x,
            prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1,
            prepk[[i]]$loc1, prepk[[i]]$n, 0, 1
          )
        } else {
          temo <- rank(-abs(Betaj), ties.method = "min")
          temo <- data.frame(temo[which(temo <= numj)], which(temo <= numj))
          temo <- temo[order(temo[, 1]), ]
          cvPL[i, ] <- cvtrimC(
            Betaj[temo[, 2]], numj, numi2, temo[, 2] - 1, prep0$x, prep0$N,
            prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x,
            prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1,
            prepk[[i]]$loc1, prepk[[i]]$n, 0, 1
          )
        }
      }
    } else {
      cvPL <- matrix(NA, nrow = nfolds, ncol = 1)
      for (i in 1:nfolds) {
        cvPL[i, ] <- cvtrimC(
          c(0.0, 0.0), 0, 0, c(0, 0), prep0$x, prep0$N, prep0$nevent,
          prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N,
          prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n,
          0, 1
        )
      }
    }

    cvraw <- cvPL / weighti
    nfoldi <- apply(!is.na(cvraw), 2, sum)
    rm(cvPL)
    cvm[[il0]] <- apply(cvraw, 2, weighted.mean, w = weighti, na.rm = TRUE)
    cv.max[il0] <- max(cvm[[il0]])

    il1 <- c(il0 - 1, il0 + 1)
    for (j in 1:2) {
      if (il1[j] >= 1 & il1[j] <= nlambdai) {
        if (is.na(cv.max[il1[j]])) {
          numi <- out$nzero[il1[j]]
          Betai <- sapply(outi, function(x) x$Beta[, il1[j]])
          Betao <- apply(Betai != 0, 2, sum)
          numi2 <- min(max(Betao), numi)

          if (numi2 > 0) {
            cvPL <- matrix(NA, nrow = nfolds, ncol = numi2)
            for (i in 1:nfolds) {
              Betaj <- Betai[, i]
              numj <- min(Betao[i], numi)
              if (numj == 0) {
                cvPL[i, ] <- cvtrimC(
                  c(0.0, 0.0), numj, numi2, c(0, 0), prep0$x, prep0$N,
                  prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n,
                  prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent,
                  prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1
                )
              } else {
                temo <- rank(-abs(Betaj), ties.method = "min")
                temo <- data.frame(temo[which(temo <= numj)], which(temo <= numj))
                temo <- temo[order(temo[, 1]), ]
                cvPL[i, ] <- cvtrimC(
                  Betaj[temo[, 2]], numj, numi2, temo[, 2] - 1, prep0$x,
                  prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n,
                  prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent,
                  prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1
                )
              }
            }
          } else {
            cvPL <- matrix(NA, nrow = nfolds, ncol = 1)
            for (i in 1:nfolds) {
              cvPL[i, ] <- cvtrimC(
                c(0.0, 0.0), 0, 0, c(0, 0), prep0$x, prep0$N,
                prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n,
                prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent,
                prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1
              )
            }
          }

          cvraw <- cvPL / weighti
          rm(cvPL)
          cvm[[il1[j]]] <- apply(cvraw, 2, weighted.mean, w = weighti,
                                 na.rm = TRUE)
          cv.max[il1[j]] <- max(cvm[[il1[j]]])
        }
      } else {
        break
      }
    }
    if (il1[j] == 1 | il1[j] == nlambdai) {
      break
    }
    if (il0 == which.max(cv.max)) {
      break
    }
    il0 <- which.max(cv.max)
  }
  index0 <- which.max(cv.max)

  Beta0 <- out$Beta[, index0]
  cuti <- which.max(cvm[[index0]])
  Beta0[abs(Beta0) <= sort(abs(Beta0), TRUE)[cuti + 1]] <- 0

  temCV0 <- data.frame(lambda = lambdai[index0], cvm = cv.max[index0], nzero = cuti)

  if (!keep.beta) {
    return(list(
      Beta = out$Beta[, indexi], Beta0 = Beta0, fit = temCV, fit0 = temCV0,
      lambda.max = lambdai[indexi], lambda.opt = lambdai[index0],
      cv.nzero = cvm[[index0]], penalty = penalty, adaptive = adaptive,
      flag = out$flag
    ))
  }
  list(
    Beta = out$Beta, Beta0 = Beta0, fit = temCV, fit0 = temCV0,
    lambda.max = lambdai[indexi], lambda.opt = lambdai[index0],
    cv.nzero = cvm[[index0]], penalty = penalty, adaptive = adaptive,
    flag = out$flag
  )
}

coxini <- function(x, y) {
  N0 <- nrow(x)
  p <- ncol(x)
  nalambda <- 10

  prep0 <- coxprep(x, y)
  wbeta <- rep(1, p)
  lambda_max <- max_lambdaC(
    prep0$x, prep0$tevent, prep0$N, prep0$nevent, prep0$nevent1,
    prep0$loc1, prep0$n, 1, wbeta, N0
  )
  lambda_min <- ifelse(N0 >= p, lambda_max * 0.0001, lambda_max * 0.01)
  alambda <- lambda_max * (lambda_min / lambda_max) ^
    (c(0:(nalambda - 1)) / (nalambda - 1))

  repeat {
    outi <- coxEnet(x, y, alpha = 0.0, lambda = alambda, isd = FALSE,
                    ifast = TRUE, keep.beta = TRUE)
    if (!is.null(outi)) {
      break
    }
    alambda <- alambda * 2.0
  }

  indexi <- ncol(outi$Beta)
  beta0 <- outi$Beta[, indexi]
  wbeta <- 1 / abs(beta0)
  sgn <- sign(beta0[1:p])
  list(wbeta = wbeta, sgn = sgn, lambda = alambda[indexi])
}

coxprep <- function(x, y) {
  N0 <- nrow(x)
  oi <- order(y[, "status"], decreasing = TRUE)
  x <- x[oi, ]
  y <- y[oi, ]
  oi <- order(y[, "time"])
  x <- x[oi, ]
  y <- y[oi, ]

  i1 <- which(y[, "status"] == 1)
  mi1 <- min(i1) - 1
  if (mi1 != 0) {
    x <- x[-c(1:mi1), ]
    y <- y[-c(1:mi1), ]
  }
  ty <- y[, "time"]
  tevent <- y[, "status"]
  N <- nrow(x)

  dty <- duplicated(ty)

  if (any(dty)) {
    tevent0 <- tevent
    tevent0[which(dty)] <- 0

    ievent <- cumsum(tevent0)
    loc1 <- which(tevent0 == 1)
    nevent <- table(ievent)
    n <- length(unique(ievent))
    nevent1 <- tapply(tevent == 1, ievent, sum)
  } else {
    ievent <- cumsum(tevent)
    loc1 <- which(tevent == 1)
    nevent <- table(ievent)
    n <- length(unique(ievent))
    nevent1 <- rep(1, n)
  }

  list(
    x = x, N0 = N0, tevent = tevent, N = N, nevent = nevent,
    nevent1 = nevent1, loc1 = loc1, n = n
  )
}

coxsplit <- function(y, nfolds) {
  N <- nrow(y)
  sample(rep(seq(nfolds), length = N))
}
