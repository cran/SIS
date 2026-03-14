loglik <- function(X, y, beta, family) {
  K <- dim(beta)[2]
  link <- cbind(1, X) %*% beta
  yrep <- repmat(y, 1, K)
  if (family == "gaussian") {
    return(apply((yrep - link)^2, 2, sum))
  }
  if (family == "poisson") {
    return(apply(exp(link) - yrep * link, 2, sum))
  }
  if (family == "binomial") {
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
  }
}

repmat <- function(X, m, n) {
  ## R equivalent of repmat (matlab)
  X <- as.matrix(X)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}

getdf <- function(coef.beta) {
  apply(abs(coef.beta) > 1e-10, 2, sum)
}


margcoef <- function(x, y, condind = NULL, family, null.model = FALSE, iterind) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  ones <- rep(1, n)
  candind <- setdiff(1:p, condind)
  if (iterind == 0) {
    if (family == "cox") {
      margcoef <- abs(cor(x, y[, 1]))
    } else if(family == 'multinom'){
      margcoef <- sapply(candind, mg, x, y, ones, family, NULL)
    } else {
      margcoef <- abs(cor(x, y))
    }
  } else {
    if (null.model == TRUE) {
      if (is.null(condind) == TRUE) {
        x <- x[sample(1:n), ]
      }
      if (is.null(condind) == FALSE) {
        x[, candind] <- x[sample(1:n), candind]
      }
    }
    margcoef <- sapply(candind, mg, x, y, ones, family, condind)
  }
  return(margcoef)
}

mg <- function(index, x = x, y = y, ones = ones, family = family, condind = condind) {
  margfit <- -switch(family,
    gaussian = glm.fit(cbind(ones, x[, index], x[, condind]), y, family = gaussian())$deviance,
    binomial = glm.fit(cbind(ones, x[, index], x[, condind]), y, family = binomial())$deviance,
    poisson = glm.fit(cbind(ones, x[, index], x[, condind]), y, family = poisson())$deviance,
    cox = -coxph(y ~ cbind(x[, index], x[, condind]))$loglik[2],
    multinom = multinom(y ~ cbind(x[, index], x[, condind]),trace = FALSE)$deviance
  )
}


obtain.ix0 <- function(x, y, s1, s2, family, nsis, iter, varISIS, perm, q, greedy, greedy.size, iterind) {
  ix0all <- NULL
  if (iter == FALSE) {
    margcoef <- margcoef(x, y, family = family, null.model = FALSE, iterind = iterind)
    rankcoef <- sort(margcoef, decreasing = TRUE, index.return = TRUE)
    ix0 <- rankcoef$ix[1:nsis]
    ix0all <- rankcoef$ix
  } else {
    if (varISIS == "vanilla") {
      margcoef <- margcoef(x, y, family = family, null.model = FALSE, iterind = iterind)
      rankcoef <- sort(margcoef, decreasing = TRUE, index.return = TRUE)
      ix0all <- rankcoef$ix
      if (perm == FALSE) {
        ix0 <- rankcoef$ix[1:floor((2 / 3) * nsis)]
      } else {
        count <- 0
        repeat {
          count <- count + 1
          randcoef <- margcoef(x, y, family = family, null.model = TRUE, iterind = iterind)
          if (length(which(margcoef >= quantile(randcoef, q))) > 0) {
            break
          }
          if (count > 10) {
            break
          }
        }
        if (greedy == FALSE) {
          if (length(which(margcoef >= quantile(randcoef, q))) >= 2) {
            length1 <- length(which(margcoef >= quantile(randcoef, q)))
            above.thresh <- rankcoef$ix[1:length1]
            ix0 <- rankcoef$ix[1:floor((2 / 3) * nsis)]
            ix0 <- sort(intersect(ix0, above.thresh))
          } else {
            ix0 <- rankcoef$ix[1:2]
          }
        } else {
          if (greedy.size == 1) {
            ix0 <- rankcoef$ix[1:2]
          } else {
            ix0 <- rankcoef$ix[1:greedy.size]
          }
        }
      }
    } else {
      if (family == "cox") {
        margcoef1 <- margcoef(x[s1, ], y[s1, ], family = family, null.model = FALSE, iterind = iterind)
        margcoef2 <- margcoef(x[s2, ], y[s2, ], family = family, null.model = FALSE, iterind = iterind)
      } else {
        margcoef1 <- margcoef(x[s1, ], y[s1], family = family, null.model = FALSE, iterind = iterind)
        margcoef2 <- margcoef(x[s2, ], y[s2], family = family, null.model = FALSE, iterind = iterind)
      }

      rankcoef1 <- sort(margcoef1, decreasing = TRUE, index.return = TRUE)
      rankcoef2 <- sort(margcoef2, decreasing = TRUE, index.return = TRUE)
      if (perm == FALSE) {
        if (varISIS == "aggr") {
          ix01 <- rankcoef1$ix[1:floor((2 / 3) * nsis)]
          ix02 <- rankcoef2$ix[1:floor((2 / 3) * nsis)]
          ix0 <- sort(intersect(ix01, ix02))
          if (length(ix0) <= 1) {
            ix0 <- int.size.k(rankcoef1$ix, rankcoef2$ix, 2)
          }
        }
        if (varISIS == "cons") {
          iensure <- intensure(floor((2 / 3) * nsis), l1 = rankcoef1$ix, l2 = rankcoef2$ix, k = floor((2 / 3) * nsis))
          ix01 <- rankcoef1$ix[1:iensure]
          ix02 <- rankcoef2$ix[1:iensure]
          ix0 <- sort(intersect(ix01, ix02))
        }
      } else {
        count <- 0
        repeat {
          count <- count + 1
          randcoef1 <- margcoef(x[s1, ], y[s1], family = family, null.model = TRUE, iterind = iterind)
          randcoef2 <- margcoef(x[s2, ], y[s2], family = family, null.model = TRUE, iterind = iterind)
          if (length(which(margcoef1 >= quantile(randcoef1, q))) > 0 && length(which(margcoef2 >= quantile(
            randcoef2,
            q
          ))) > 0) {
            break
          }
          if (count > 10) break
        }
        if (greedy == FALSE) {
          length1 <- length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 <- length(which(margcoef2 >= quantile(randcoef2, q)))
          above.thresh.1 <- rankcoef1$ix[1:length1]
          above.thresh.2 <- rankcoef2$ix[1:length2]
          ix01 <- rankcoef1$ix[1:floor((2 / 3) * nsis)]
          ix02 <- rankcoef2$ix[1:floor((2 / 3) * nsis)]
          ix01 <- sort(intersect(ix01, above.thresh.1))
          ix02 <- sort(intersect(ix02, above.thresh.2))
          ix0 <- sort(intersect(ix01, ix02))
          if (length(ix0) <= 1) {
            ix0 <- int.size.k(rankcoef1$ix, rankcoef2$ix, 2)
          }
        } else {
          if (greedy.size == 1) {
            ix0 <- int.size.k(rankcoef1$ix, rankcoef2$ix, 2)
          } else {
            ix0 <- int.size.k(rankcoef1$ix, rankcoef2$ix, greedy.size)
          }
        }
      }
    }
  }
  return(list(ix0 = ix0, ix0all = ix0all))
}

obtain.newix <- function(x, y, ix1, candind, s1, s2, family, pleft, varISIS, perm, q, greedy, greedy.size, iterind) {
  newixall <- NULL
  if (varISIS == "vanilla") {
    margcoef <- margcoef(x, y, ix1, family = family, null.model = FALSE, iterind = iterind)
    rankcoef <- sort(margcoef, decreasing = TRUE, index.return = TRUE)
    newixall <- candind[rankcoef$ix]
    if (perm == FALSE) {
      if (pleft > 0) {
        newix <- candind[rankcoef$ix[1:pleft]]
      } else {
        newix <- NULL
      }
    } else {
      randcoef <- margcoef(x, y, ix1, family = family, null.model = TRUE, iterind = iterind)
      if (length(which(margcoef >= quantile(randcoef, q))) > 0) {
        if (greedy == FALSE) {
          length1 <- length(which(margcoef >= quantile(randcoef, q)))
          above.thresh <- candind[rankcoef$ix[1:length1]]
          newix <- candind[rankcoef$ix[1:pleft]]
          newix <- sort(intersect(newix, above.thresh))
        } else {
          newix <- candind[rankcoef$ix[1:greedy.size]]
        }
      } else {
        newix <- NULL
      }
    }
  } else {
    margcoef1 <- margcoef(x[s1, ], y[s1], ix1, family = family, null.model = FALSE, iterind = iterind)
    margcoef2 <- margcoef(x[s2, ], y[s2], ix1, family = family, null.model = FALSE, iterind = iterind)
    rankcoef1 <- sort(margcoef1, decreasing = TRUE, index.return = TRUE)
    rankcoef2 <- sort(margcoef2, decreasing = TRUE, index.return = TRUE)
    if (perm == FALSE) {
      if (pleft > 0) {
        if (varISIS == "aggr") {
          newix1 <- candind[rankcoef1$ix[1:pleft]]
          newix2 <- candind[rankcoef2$ix[1:pleft]]
          newix <- sort(intersect(newix1, newix2))
        }
        if (varISIS == "cons") {
          iensure <- intensure(pleft, l1 = rankcoef1$ix, l2 = rankcoef2$ix, k = pleft)
          newix1 <- candind[rankcoef1$ix[1:iensure]]
          newix2 <- candind[rankcoef2$ix[1:iensure]]
          newix <- sort(intersect(newix1, newix2))
        }
      } else {
        newix <- NULL
      }
    } else {
      randcoef1 <- margcoef(x[s1, ], y[s1], ix1, family = family, null.model = TRUE, iterind = iterind)
      randcoef2 <- margcoef(x[s2, ], y[s2], ix1, family = family, null.model = TRUE, iterind = iterind)
      if (length(which(margcoef1 >= quantile(randcoef1, q))) > 0 && length(which(margcoef2 >= quantile(randcoef2, q))) > 0) {
        if (greedy == FALSE) {
          length1 <- length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 <- length(which(margcoef2 >= quantile(randcoef2, q)))
          above.thresh.1 <- candind[rankcoef1$ix[1:length1]]
          above.thresh.2 <- candind[rankcoef2$ix[1:length2]]
          newix1 <- candind[rankcoef1$ix[1:pleft]]
          newix2 <- candind[rankcoef2$ix[1:pleft]]
          newix1 <- sort(intersect(newix1, above.thresh.1))
          newix2 <- sort(intersect(newix2, above.thresh.2))
          newix <- sort(intersect(newix1, newix2))
        } else {
          length1 <- length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 <- length(which(margcoef2 >= quantile(randcoef2, q)))
          newix1 <- candind[rankcoef1$ix[1:length1]]
          newix2 <- candind[rankcoef2$ix[1:length2]]
          iensure <- intensure(greedy.size, l1 = newix1, l2 = newix2, k = greedy.size)
          if (is.null(iensure)) {
            newix <- NULL
          } else {
            newix <- sort(intersect(newix1[1:iensure], newix2[1:iensure]))
          }
        }
      } else {
        newix <- NULL
      }
    }
  }
  return(list(newix = newix, newixall = newixall))
}


intensure <- function(i, l1, l2, k) {
  for (j in i:length(l1)) {
    if (length(intersect(l1[1:j], l2[1:j])) >= k) {
      return(j)
    }
  }
  # if (length(intersect(l1[1:i], l2[1:i])) >= k)
  #    return(i) else return(intensure(i + 1, l1, l2, k))
}

int.size.k <- function(l1, l2, k) {
  iensure <- intensure(k, l1 = l1, l2 = l2, k = k)
  ix01 <- l1[1:iensure]
  ix02 <- l2[1:iensure]
  ix0 <- sort(intersect(ix01, ix02))
  return(ix0)
}

calculate.nsis <- function(family, varISIS, n, p) {
  if (varISIS == "aggr") {
    nsis <- floor(n / log(n))
  } else {
    if (family == "gaussian") {
      nsis <- floor(n / log(n))
    }
    if (family == "binomial" | family == 'multinom') {
      nsis <- floor(n / (4 * log(n)))
    }
    if (family == "poisson") {
      nsis <- floor(n / (2 * log(n)))
    }
    if (family == "cox") {
      nsis <- floor(n / (4 * log(n)))
    }
  }
  if (p < n) {
    nsis <- p
  }
  return(nsis)
}
