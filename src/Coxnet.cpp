// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X) {
  Eigen::VectorXd mX = X.colwise().mean(), sdX(X.cols()), sdXi;
  X.rowwise() -= mX.transpose();
  sdX = X.colwise().norm() / sqrt((double)X.rows());
  sdXi = 1.0 / sdX.array();
  X = X * sdXi.asDiagonal();
  return List::create(Named("x") = X, Named("sd") = sdX);
}

/*****  Lambda path (max)  *****/
// [[Rcpp::export]]
double max_lambdaC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N,
                   Eigen::VectorXi nevent, Eigen::VectorXi nevent1,
                   Eigen::VectorXi loc1, int n, double alpha,
                   Eigen::VectorXd wbeta, int N0) {
  int i, j, q;
  double denS = N, c1 = 0.0;
  Eigen::VectorXd Li(N), lli(N);

  for (i = 0; i < n; i++) {
    c1 += (nevent1(i) / denS);
    denS -= nevent(i);
    for (j = loc1(i) - 1, q = 0; q < nevent(i); j++, q++) {
      lli(j) = tevent(j) - c1;
    }
  }
  Li = (lli.transpose() * X).cwiseAbs() / N0;
  Li = Li.array() / wbeta.array() / alpha;
  return Li.maxCoeff();
}

/*****  Derivatives of log-pl of eta (1st&2nd order), ties  *****/
void dletaCm(Eigen::VectorXd & exb, Eigen::VectorXd & tevent, int & N,
             Eigen::VectorXi & nevent, Eigen::VectorXi & nevent1,
             Eigen::VectorXi & loc1, int & n, Eigen::VectorXd & pl1,
             Eigen::VectorXd & pl2, int & ifast, int & itwo) {
  int i, j, q, ipl2 = 0;
  double denSi, c1 = 0.0, c2 = 0.0;
  Eigen::VectorXd denS(n);

  if (ifast == 0 || itwo == 1) goto two;
  denSi = exb.sum();
  for (i = 0; i < n; ++i) {
    c1 += (nevent1(i) / denSi);
    c2 += (nevent1(i) / pow(denSi, 2));
    for (j = loc1(i) - 1, q = 0; q < nevent(i); j++, q++) {
      denSi -= exb(j);
      pl1(j) = tevent(j) - exb(j) * c1;
      pl2(j) = exb(j) * (c1 - exb(j) * c2);
      if (pl2(j) <= 0.0) ipl2 = 1;
    }
  }
  if (ipl2 == 1) {
    itwo = 1;
    if (ifast == 0) {
      goto two;
    }
  }
  return;

  two:
  denSi = 0.0;
  c1 = 0.0;
  c2 = 0.0;
  for (i = n - 1; i >= 0; --i) {
    for (j = loc1(i) - 1, q = 0; q < nevent(i); j++, q++) {
      denSi += exb(j);
    }
    denS(i) = denSi;
  }
  for (i = 0; i < n; ++i) {
    c1 += (nevent1(i) / denS(i));
    c2 += (nevent1(i) / pow(denS(i), 2));
    for (j = loc1(i) - 1, q = 0; q < nevent(i); j++, q++) {
      pl1(j) = tevent(j) - exb(j) * c1;
      pl2(j) = exb(j) * (c1 - exb(j) * c2);
    }
  }
}

/*****  Log-pl of eta, ties  *****/
double pletaCm(Eigen::VectorXd & xb, Eigen::VectorXd & exb,
               Eigen::VectorXi & nevent, Eigen::VectorXi & nevent1,
               Eigen::VectorXi & loc1, int & n, int & ifast, int & itwo) {
  int i, j, q, iSS = 0;
  double ll = 0.0, SSi;
  Eigen::VectorXd SS(n);

  if (ifast == 0 || itwo == 1) goto two;
  SSi = exb.sum();
  for (i = 0; i < n; ++i) {
    if (SSi <= 0.0) iSS = 1;
    for (j = loc1(i) - 1, q = 0; q < nevent1(i); j++, q++) {
      ll += xb(j);
    }
    ll -= nevent1(i) * log(SSi);
    for (j = loc1(i) - 1, q = 0; q < nevent(i); j++, q++) {
      SSi -= exb(j);
    }
  }
  if (iSS == 1) {
    itwo = 1;
    if (ifast == 0) {
      goto two;
    }
  }
  return ll;

  two:
  ll = 0.0;
  SSi = 0.0;
  for (i = n - 1; i >= 0; --i) {
    for (j = loc1(i) - 1, q = 0; q < nevent(i); j++, q++) {
      SSi += exb(j);
    }
    SS(i) = SSi;
  }
  for (i = 0; i < n; ++i) {
    for (j = loc1(i) - 1, q = 0; q < nevent1(i); j++, q++) {
      ll += xb(j) - log(SS(i));
    }
  }
  return ll;
}

/*****  Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvtrimC(Eigen::VectorXd beta, int nn, int nn2,
                        Eigen::VectorXi loco, Eigen::MatrixXd XF, int NF,
                        Eigen::VectorXi neventF, Eigen::VectorXi nevent1F,
                        Eigen::VectorXi loc1F, int nF, Eigen::MatrixXd X, int N,
                        Eigen::VectorXi nevent, Eigen::VectorXi nevent1,
                        Eigen::VectorXi loc1, int n, int ifast, int itwo) {
  int i, j;
  double lli, lfi;
  Eigen::VectorXd cv, xb = Eigen::VectorXd::Zero(N),
                  xbF = Eigen::VectorXd::Zero(NF);
  Eigen::VectorXd exb(N), exbF(NF);

  if (nn2 > 0) {
    cv = Eigen::VectorXd::Zero(nn2);

    if (nn == 0) {
      exb = (xb.array()).exp();
      lli = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      exbF = (xbF.array()).exp();
      lfi = pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
      cv(0) = lfi - lli;
    } else {
      for (i = 0; i < nn; i++) {
        j = loco(i);
        xb += X.col(j) * beta(i);
        exb = (xb.array()).exp();
        lli = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        xbF += XF.col(j) * beta(i);
        exbF = (xbF.array()).exp();
        lfi = pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
        cv(i) = lfi - lli;
      }
    }
    if (nn2 > nn) {
      for (i = nn; i < nn2; i++) {
        cv(i) = cv(nn - 1);
      }
    }
  } else {
    cv = Eigen::VectorXd::Zero(1);

    exb = (xb.array()).exp();
    lli = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
    exbF = (xbF.array()).exp();
    lfi = pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
    cv(0) = lfi - lli;
  }

  return cv;
}

/*****  Enet (L1+L2)  *****/
// [[Rcpp::export]]
List coxenetC(Eigen::MatrixXd X, Eigen::VectorXd tevent, double alpha,
              Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
              int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1,
              Eigen::VectorXi loc1, int n, int p, int N0, double thresh,
              int maxit, int ifast) {
  int i, j, it, il, iadd, ia = 0, itwo = 0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2,
    objQi, objQj;
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active = Eigen::VectorXi::Zero(p),
                  iactive = Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag = Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb = Eigen::VectorXd::Constant(N, 1.0),
                  xb = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);

  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0 = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0 = -ll0 / N0;

  for (il = 0; il < nlambda; ++il) {
    lambda1 = lambda(il) * alpha * wbeta;
    lambda2 = lambda(il) * (1 - alpha);
    lambda1i = N0 * lambda1;
    lambda2i = N0 * lambda2;

    for (i = 0; i < p; ++i) {
      if (iactive(i) == 0) {
        PLi = pl1.dot(X.col(i));
        if (std::abs(PLi) > lambda1i(i)) {
          active(ia) = i;
          iactive(i) = 1;
          ++ia;
        }
      }
    }

    it = 0;
    local:
    while (1) {
      ++it;

      objQi = 0.0;
      objQj = 0.0;
      for (i = 0; i < ia; ++i) {
        j = active(i);
        PLi2 = pl2.dot(X.col(j).cwiseAbs2());
        zi = beta0(j) * PLi2 + pl1.dot(X.col(j));
        if (zi > lambda1i(j)) {
          b0 = (zi - lambda1i(j)) / (lambda2i + PLi2);
          db0 = beta0(j) - b0;
          beta0(j) = b0;
          pl1 += (pl2.cwiseProduct(X.col(j))) * db0;
          xb -= db0 * X.col(j);
          objQj += std::abs(b0) * lambda1(j);
          objQi += pow(b0, 2);
        } else if (zi < -lambda1i(j)) {
          b0 = (zi + lambda1i(j)) / (lambda2i + PLi2);
          db0 = beta0(j) - b0;
          beta0(j) = b0;
          pl1 += (pl2.cwiseProduct(X.col(j))) * db0;
          xb -= db0 * X.col(j);
          objQj += std::abs(b0) * lambda1(j);
          objQi += pow(b0, 2);
        } else {
          b0 = 0.0;
          if (beta0(j) != b0) {
            db0 = beta0(j) - b0;
            beta0(j) = b0;
            pl1 += (pl2.cwiseProduct(X.col(j))) * db0;
            xb -= db0 * X.col(j);
          }
        }
      }

      ll1 = ll0;
      obj1 = obj0;
      exb = (xb.array()).exp();
      ll0 = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if (ifast == 1 && itwo == 1) goto exit;
      obj0 = -ll0 / N0 + objQj + objQi * lambda2 / 2.0;

      if (std::abs(ll1 - ll0) < std::abs(thresh * ll1)) {
        flag(il) = 0;
        break;
      }
      if (std::abs(obj1 - obj0) < std::abs(thresh * obj1)) {
        flag(il) = 0;
        break;
      }
      if (obj0 != obj0) {
        flag(il) = 2;
        goto exit;
      }
      if (it >= maxit) {
        flag(il) = 1;
        goto exit;
      }

      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if (ifast == 1 && itwo == 1) goto exit;
    }

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if (ifast == 1 && itwo == 1) goto exit;
    iadd = 0;
    for (i = 0; i < p; ++i) {
      if (iactive(i) == 0) {
        PLi = pl1.dot(X.col(i));
        if (std::abs(PLi) > lambda1i(i)) {
          active(ia) = i;
          iactive(i) = 1;
          ++ia;
          iadd = 1;
        }
      }
    }
    if (iadd == 1) goto local;

    locbeta(il) = ll0;
    Betas.col(il) = beta0;
  }

  exit:
  if (ifast == 1 && itwo == 1 && il > 0) --il;
  return List::create(
    Named("Beta") = Betas, Named("flag") = flag, Named("ll") = locbeta,
    Named("nlambda") = il
  );
}

/*****  Enet (L1+L2)  *****/
/*****  cross-validation PL  *****/
// [[Rcpp::export]]
List cvcoxenetC(Eigen::MatrixXd X, Eigen::VectorXd tevent, double alpha,
                Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
                int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1,
                Eigen::VectorXi loc1, int n, int p, int N0, double thresh,
                int maxit, int ifast, Eigen::MatrixXd XF, int NF,
                Eigen::VectorXi neventF, Eigen::VectorXi nevent1F,
                Eigen::VectorXi loc1F, int nF) {
  int i, j, it, il, iadd, ia = 0, itwo = 0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2,
    objQi, objQj;
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active = Eigen::VectorXi::Zero(p),
                  iactive = Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag = Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb = Eigen::VectorXd::Constant(N, 1.0),
                  xb = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);

  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0 = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0 = -ll0 / N0;

  for (il = 0; il < nlambda; ++il) {
    lambda1 = lambda(il) * alpha * wbeta;
    lambda2 = lambda(il) * (1 - alpha);
    lambda1i = N0 * lambda1;
    lambda2i = N0 * lambda2;

    for (i = 0; i < p; ++i) {
      if (iactive(i) == 0) {
        PLi = pl1.dot(X.col(i));
        if (std::abs(PLi) > lambda1i(i)) {
          active(ia) = i;
          iactive(i) = 1;
          ++ia;
        }
      }
    }

    it = 0;
    local:
    while (1) {
      ++it;

      objQi = 0.0;
      objQj = 0.0;
      for (i = 0; i < ia; ++i) {
        j = active(i);
        PLi2 = pl2.dot(X.col(j).cwiseAbs2());
        zi = beta0(j) * PLi2 + pl1.dot(X.col(j));
        if (zi > lambda1i(j)) {
          b0 = (zi - lambda1i(j)) / (lambda2i + PLi2);
          db0 = beta0(j) - b0;
          beta0(j) = b0;
          pl1 += (pl2.cwiseProduct(X.col(j))) * db0;
          xb -= db0 * X.col(j);
          objQj += std::abs(b0) * lambda1(j);
          objQi += pow(b0, 2);
        } else if (zi < -lambda1i(j)) {
          b0 = (zi + lambda1i(j)) / (lambda2i + PLi2);
          db0 = beta0(j) - b0;
          beta0(j) = b0;
          pl1 += (pl2.cwiseProduct(X.col(j))) * db0;
          xb -= db0 * X.col(j);
          objQj += std::abs(b0) * lambda1(j);
          objQi += pow(b0, 2);
        } else {
          b0 = 0.0;
          if (beta0(j) != b0) {
            db0 = beta0(j) - b0;
            beta0(j) = b0;
            pl1 += (pl2.cwiseProduct(X.col(j))) * db0;
            xb -= db0 * X.col(j);
          }
        }
      }

      ll1 = ll0;
      obj1 = obj0;
      exb = (xb.array()).exp();
      ll0 = pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if (ifast == 1 && itwo == 1) goto exit;
      obj0 = -ll0 / N0 + objQj + objQi * lambda2 / 2.0;

      if (std::abs(ll1 - ll0) < std::abs(thresh * ll1)) {
        flag(il) = 0;
        break;
      }
      if (std::abs(obj1 - obj0) < std::abs(thresh * obj1)) {
        flag(il) = 0;
        break;
      }
      if (obj0 != obj0) {
        flag(il) = 2;
        goto exit;
      }
      if (it >= maxit) {
        flag(il) = 1;
        goto exit;
      }

      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if (ifast == 1 && itwo == 1) goto exit;
    }

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if (ifast == 1 && itwo == 1) goto exit;
    iadd = 0;
    for (i = 0; i < p; ++i) {
      if (iactive(i) == 0) {
        PLi = pl1.dot(X.col(i));
        if (std::abs(PLi) > lambda1i(i)) {
          active(ia) = i;
          iactive(i) = 1;
          ++ia;
          iadd = 1;
        }
      }
    }
    if (iadd == 1) goto local;

    locbeta(il) = ll0;
    Betas.col(il) = beta0;

    xbF = Eigen::VectorXd::Zero(NF);
    for (i = 0; i < ia; i++) {
      j = active(i);
      xbF += XF.col(j) * beta0(j);
    }
    exbF = (xbF.array()).exp();
    locbetaF(il) = pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }

  exit:
  if (ifast == 1 && itwo == 1 && il > 0) --il;
  return List::create(
    Named("Beta") = Betas, Named("flag") = flag, Named("ll") = locbeta,
    Named("lf") = locbetaF, Named("nlambda") = il
  );
}
