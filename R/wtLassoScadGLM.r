INDEPtuneSCADglm <- function (x, y, wt.initsoln = wt.initsoln, xtune, ytune, dlambda = NULL,
    maxlambda=NULL, family = family, nopenalty.subset = NULL) 
{
   if(is.null(maxlambda)) maxlambda = findSCADMaxLambda(x, y, wt.initsoln = wt.initsoln, family = family)
        
   lambda.cand = maxlambda * (2^(seq(0, -20, -0.2)))

   oldsoln = NULL
   tuneer = NULL
   
   x=as.matrix(x)
   weight=rep(1, nrow(x))
   offset=rep(0, nrow(x))
   
   iloop = 1    
   oldsoln[[iloop]] = scadglm(x, y, wt.initsoln = wt.initsoln, 
        lambda = lambda.cand[iloop], initsoln = NULL, family = family, 
        weight = weight, offset = offset, nopenalty.subset = nopenalty.subset)
   tuneer[iloop] =familyglmNEGloglik(xtune, ytune, oldsoln[[iloop]], family = family)

   for (iloop in 2:length(lambda.cand)) {
       oldsoln[[iloop]] = scadglm(x, y, wt.initsoln = wt.initsoln, 
           lambda = lambda.cand[iloop], initsoln = oldsoln[[iloop - 
           1]], family = family, weight = weight, offset = offset, 
           nopenalty.subset = nopenalty.subset)
       tuneer[iloop] = familyglmNEGloglik(xtune, ytune, 
           oldsoln[[iloop]], family = family)
    }
    best.lambda.ind = which.min(tuneer)
    w = scadglm(x, y, wt.initsoln = wt.initsoln, lambda = lambda.cand[best.lambda.ind], 
        initsoln = oldsoln[[best.lambda.ind]], family = family, 
        nopenalty.subset = nopenalty.subset)
   return(list(wt.initsoln = wt.initsoln, tuneer = tuneer, best.lambda.ind = best.lambda.ind, 
        maxlambda = maxlambda, w = w, used.lambda = lambda.cand[1:iloop], 
        best.lambda = lambda.cand[best.lambda.ind]))
}


############################################################################################################################################################
positivepart <- function(x) {
  return(ifelse(x>=0, x, 0))
}
############################################################################################################################################################
scadderiv <- function(ftheta, fa, flambda) {
 return(flambda*(1-(1-apply(as.matrix(fa*flambda-abs(ftheta)), 1, positivepart)/((fa-1)*flambda))*as.numeric(abs(ftheta)>flambda)))
}
############################################################################################################################################################
familyglmNEGloglik <-function(x, y, fwcoef, family=family, weight=NULL, offset=NULL) {
  # this function returns the loglikelihood for GLMs
  if(is.null(offset)) offset=rep(0, nrow(x))
  if(is.null(weight)) weight=rep(1, nrow(x))
   eta=as.matrix(x)%*%as.matrix(fwcoef[-1])+fwcoef[1]+offset
   return(switch(family$family,
        gaussian = sum(((y-eta)^2)*weight),
        binomial = sum(-(y*eta-log(1+exp(eta)))*weight),
        poisson = sum(-(y*eta-exp(eta))*weight)))
}
############################################################################################################################################################
findSCADMaxLambda <- function(x, y, wt.initsoln=wt.initsoln, family=family) {
x=cbind(1, x)
n=nrow(x)
p=ncol(x)-1

ones=rep(1, n)
intercept.coef=coef(glm.fit(ones, y, intercept=FALSE, family=family))

variance <- family$variance
linkinv <- family$linkinv
mu.eta <- family$mu.eta
dev.resids <- family$dev.resids

eta <- intercept.coef*ones  # +offset
mu <- linkinv(eta)
mu.eta.eta <- mu.eta(eta)
w <- (mu.eta.eta^2/variance(mu))^0.5
z <- (y-mu)/mu.eta.eta
xw <- x * w
wz <- w * z
corr <-  drop(t(wz) %*% xw)
temp.all.lambda=apply(cbind(abs(corr)/n, abs(wt.initsoln)), 1, scadafindmaxlam)
lambda=max(temp.all.lambda[-1])
return(lambda)
}
############################################################################################################################################################
"scadafindmaxlam" <-function(temp) {
 return(findmaxlam(temp, 3.7))
}
############################################################################################################################################################
"findmaxlam" <- function(wux, a) {
 return(ifelse(wux[1]>=wux[2], wux[1], wux[1]+(wux[2]-wux[1])/a))
}
############################################################################################################################################################
CVscadglm <- function(x, y, wt.initsoln=NULL, dlambda=NULL, maxlambda=NULL, folds=NULL, 
  family=family, weight=NULL, offset=NULL, nopenalty.subset=NULL) {

x=as.matrix(x)
n=nrow(x)
if(is.null(offset)) offset=rep(0, nrow(x))
if(is.null(weight)) weight=rep(1, nrow(x))
if(is.null(wt.initsoln)) wt.initsoln=rep(0, ncol(x)+1)
if(is.null(maxlambda)) maxlambda=findSCADMaxLambda(x, y, wt.initsoln=wt.initsoln, family=family)

  lambda.cand=maxlambda*(2^(seq(0, -20, -0.2)))
 
if(is.null(folds)) {
temp= sample(1:n, n, replace = FALSE)
kfold=10
for(i in 1:kfold){
  folds[[i]]=setdiff(1:n, temp[seq(i, n, kfold)])
  }
} 
 
n.folds=length(folds)
nobs.fold=NULL
for (i in 1:n.folds) nobs.fold[i]=length(folds[[i]])

oldsoln=NULL
oldsoln.save=NULL

tuneer=NULL
iloop=1
tuneer[iloop]=0
for (j in 1:n.folds) {
  oldsoln[[j]]=scadglm(x[folds[[j]],], y[folds[[j]]], wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=NULL, family = family, weight = weight[folds[[j]]], offset = offset[folds[[j]]], nopenalty.subset=nopenalty.subset)
  tuneer[iloop]=tuneer[[iloop]]+familyglmNEGloglik(x[-folds[[j]],], y[-folds[[j]]], oldsoln[[j]], family=family, weight=weight[-folds[[j]]], offset=offset[-folds[[j]]])
}
oldsoln.save[[iloop]]=oldsoln[[1]]

temp.initsoln=wt.initsoln[-1]

for (iloop in 2:length(lambda.cand)) {
  if((max(abs(temp.initsoln))>0) && (min(abs(temp.initsoln[abs(temp.initsoln)>0]))>3.7*lambda.cand[iloop-1])) break
  tuneer[iloop]=0
  for (j in 1:n.folds) {
    oldsoln[[j]]=scadglm(x[folds[[j]],], y[folds[[j]]], wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=oldsoln[[j]], family = family, weight = weight[folds[[j]]], offset = offset[folds[[j]]], nopenalty.subset=nopenalty.subset)
    tuneer[iloop]=tuneer[[iloop]]+familyglmNEGloglik(x[-folds[[j]],], y[-folds[[j]]], oldsoln[[j]], family=family, weight=weight[-folds[[j]]], offset=offset[-folds[[j]]])
  }
oldsoln.save[[iloop]]=oldsoln[[1]]
}
best.lambda.ind=which.min(tuneer) 
w=scadglm(x, y, wt.initsoln=wt.initsoln, lambda=lambda.cand[best.lambda.ind], initsoln=oldsoln.save[[best.lambda.ind]], family = family, weight = weight, offset = offset, nopenalty.subset=nopenalty.subset)
return(list(wt.initsoln=wt.initsoln, tuneer=tuneer, best.lambda.ind=best.lambda.ind, maxlambda=maxlambda, w=w, used.lambda=lambda.cand[1:iloop], best.lambda=lambda.cand[best.lambda.ind]))
}    # end of CVscadglm function
############################################################################################################################################################
AICBICscadglm <- function(x, y, wt.initsoln=wt.initsoln, dlambda=NULL, maxlambda=NULL, 
   family=family, AICBIC='AIC', weight=NULL, offset=NULL, nopenalty.subset=NULL,  eps0=1e-3) {

x=as.matrix(x)
if(is.null(offset)) offset=rep(0, nrow(x))
if(is.null(weight)) weight=rep(1, nrow(x))
if(is.null(maxlambda)) maxlambda=findSCADMaxLambda(x, y, wt.initsoln=wt.initsoln, family=family)


  lambda.cand=maxlambda*(2^(seq(0, -20, -0.2)))
 
AB.factor=switch(AICBIC, AIC=2, BIC=log(nrow(x)))

oldsoln=NULL

tuneer=NULL
iloop=1
oldsoln[[iloop]]=scadglm(x, y, wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=NULL, family = family, weight = weight, offset = offset, nopenalty.subset=nopenalty.subset)
tuneer[iloop]=length(which(abs(oldsoln[[iloop]][-1])>eps0))*AB.factor+2*familyglmNEGloglik(x, y, oldsoln[[iloop]], family=family, weight=weight, offset=offset)

for (iloop in 2:length(lambda.cand)) {
 oldsoln[[iloop]]=scadglm(x, y, wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=oldsoln[[iloop-1]], family = family, weight = weight, offset = offset, nopenalty.subset=nopenalty.subset)
 tuneer[iloop]=length(which(abs(oldsoln[[iloop]][-1])>eps0))*AB.factor+2*familyglmNEGloglik(x, y, oldsoln[[iloop]], family=family, weight=weight, offset=offset)
}

best.lambda.ind=which.min(tuneer) 
w=scadglm(x, y, wt.initsoln=wt.initsoln, lambda=lambda.cand[best.lambda.ind], initsoln=oldsoln[[best.lambda.ind]], family = family, weight = weight, offset = offset, nopenalty.subset=nopenalty.subset)
return(list(wt.initsoln=wt.initsoln, tuneer=tuneer, best.lambda.ind=best.lambda.ind, maxlambda=maxlambda, w=w, used.lambda=lambda.cand[1:iloop], best.lambda=lambda.cand[best.lambda.ind]))
}    # end of CVscadglm function
############################################################################################################################################################
"fullscadglm" <- function(x, y, lambda, initsoln=NULL, family = binomial(), 
   weight = NULL, offset = NULL, function.precision=1e-8, nopenalty.subset=NULL, eps0=1e-6) {

  if(is.null(offset)) offset=rep(0, nrow(x))
  if(is.null(weight)) weight=rep(1, nrow(x))
  
nobs=nrow(x)
pdim=ncol(x)
if(is.null(initsoln)) initsoln=rep(0, pdim+1)

wold=initsoln
test=1
while(test) {
lassoweight=nobs*scadderiv(abs(wold[-1]), 3.7, lambda)
lassoweight[nopenalty.subset]=0
wnew=wtlassoglm(x, y, lassoweight=lassoweight, initsoln=wold, family=family, weight=weight, 
      offset=offset, lambda2=0, function.precision=function.precision)$w
if(sum((wold-wnew)^2)<eps0) test=0
 wold=wnew
}
      
return(wnew)
}
############################################################################################################################################################
"scadglm" <- function(x, y, wt.initsoln=NULL, lambda, initsoln=NULL, family = binomial(), 
   weight = NULL, offset = NULL, function.precision=1e-8, nopenalty.subset=NULL) {

   x=as.matrix(x)
  if(is.null(offset)) offset=rep(0, nrow(x))
  if(is.null(weight)) weight=rep(1, nrow(x))
  
nobs=nrow(x)
pdim=ncol(x)
if(is.null(wt.initsoln)) wt.initsoln=rep(0,pdim+1)
lassoweight=nobs*scadderiv(abs(wt.initsoln[-1]), 3.7, lambda)
lassoweight[nopenalty.subset]=0
return(wtlassoglm(x, y, lassoweight=lassoweight, initsoln=initsoln, family=family, weight=weight, offset=offset, lambda2=0, function.precision=function.precision)$w)
}
############################################################################################################################################################
"wtlassoglm" <- function(x, y, lassoweight=NULL, initsoln=NULL,   family = binomial(), weight = NULL, 
      offset = NULL, lambda2=0, function.precision=1e-8) {

  if(is.null(offset)) offset=rep(0, nrow(x))
  if(is.null(weight)) weight=rep(1, nrow(x))
  
call <- match.call()
if (is.character(family))
  family <- get(family, mode = "function", envir = parent.frame())
if (is.function(family))
   family <- family()
else if (family$family=="gaussian") family <- gaussian()
else if (family$family=="binomial") family <- binomial()
else if (family$family=="poisson") family <- poisson()
if (is.null(family$family)) {
   print(family)
   stop("'family' not recognized")
}
nobs <- nrow(x)
m <- ncol(x)

if(is.null(lassoweight)) lassoweight=rep(1/nobs, m)
if(is.null(initsoln))  initsoln=rep(0, m+1)
if (family$family=="binomial") {
   if (any(y==-1)) {
       y[y==-1] <- 0
       cat("y=-1 converted to y=0.\n")
   }
}
if (length(weight) != nobs) {
   stop("Length of the weight vector != sample size")
}
if (any(weight < 0)) {
   stop("Negative weights are not allowed.")
}
if (length(offset) != nobs) {
   stop("Length of the offset vector != sample size")
}

x=cbind(1, x)

if (is.null(dimnames(x)[[2]])) {
   xnames <- c("Intercept", paste("x",seq(m),sep=""))
}
else xnames <- c("Intercept", dimnames(x)[[2]][-1])


lambda=0  

b0=initsoln
a0=abs(b0[-1])

tmpa=1:ncol(x)
force.active=1

    p <- length(tmpa)
    fk <- length(force.active)
    k <- p - fk
      param <- c(b0[tmpa], a0, 0, -b0[tmpa[-force.active]]-a0, b0[tmpa[-force.active]]-a0)
      xa <- x[ ,tmpa,drop=FALSE]
      xstate <- rep(2, p+3*k+1)
      xstate[param==0] <- 0
      dstr <- switch(family$family, gaussian=0, binomial=1, poisson=2)
      lenz <- 10+(p+3)*nobs+k
      zsmall <- rep(0, lenz)
      zsmall[1:6] <- c(nobs, lambda, lambda2/nobs, dstr, k, function.precision)
      zsmall[10 + seq((p+3)*nobs)] <- c(as.vector(xa), y, weight, offset)
      zsmall[10+(p+3)*nobs+seq(k)] = lassoweight/nobs 
      
         ptm <- proc.time()
      sol <- .Fortran("solutionwu",
                      k = as.integer(k),
                      n = as.integer(p+k),
                      nb = as.integer(p+3*k+1),
                      ne = as.integer(p+3*k),
                      hs = as.integer(xstate),
                      xn = as.double(param),
                      zsmall = as.double(zsmall),
                      lenz = as.integer(lenz),
                      inform = as.integer(0))
                      proc.time() - ptm
      b0[tmpa] <- sol$xn[1:p]
      if (k > 0) a0 <- sol$xn[(p+1):(p+k)]


      if (sol$inform != 0) warning("\nconvergence warning in corrector step\n")



    object <- list(call = call, lambda2=lambda2,  xnames = xnames, family = family, weight = weight, offset = offset, lassoweight=lassoweight, initsoln=initsoln, w=b0)

    object


    }
############################################################################################################################################################


################################################################
INDEPgetfinalSCADcoef <- function(x, y, pickind, folds=NULL,  xtune, ytune, family=binomial(),  inittype='NoPen', eps0=1e-3, detailed=FALSE) {

fp=ncol(x)
fn=nrow(x)
ones=rep(1, fn)


if(inittype=='L1') {
  L1result=CVscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, folds=folds, family=family)
  wt.initsoln=L1result$w
}
if(inittype=='NoPen') {
wt.initsoln=coef(glm.fit(cbind(ones, x[, pickind]), y, family=family))
}

SCADresult=INDEPtuneSCADglm(x[, pickind], y, wt.initsoln=wt.initsoln, xtune[, pickind], ytune, family=family)
 
tempw=SCADresult$w
tempw[abs(tempw)<eps0]=0
SCADcoef=rep(0, fp+1)
SCADcoef[c(1, pickind+1)]=tempw # SCADresult$w
  
 if (detailed) {
        return(list(wt.initsoln = wt.initsoln, SCADcoef = SCADcoef, SCADresult=SCADresult))
    }
    else {
        return(list(wt.initsoln = wt.initsoln, SCADcoef = SCADcoef))
    }  
}



######################################################################################################################

getfinalSCADcoef <- function(x, y, pickind, folds=NULL, eps0=1e-3, family=binomial(), tune.method="AIC", inittype='NoPen', detailed=FALSE) {

fp=ncol(x)
fn=nrow(x)
ones=rep(1, fn)


if(inittype=='L1') {
  L1result=CVscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, folds=folds, family=family)
  wt.initsoln=L1result$w
}
if(inittype=='NoPen') {
wt.initsoln=coef(glm.fit(cbind(ones, x[, pickind]), y, family=family))
}

SCADresult=switch(tune.method,
CV=CVscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, folds=folds, family=family),
AIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC='AIC', eps0=eps0),
BIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC='BIC', eps0=eps0))
         

tempw=SCADresult$w
tempw[abs(tempw)<eps0]=0
SCADcoef=rep(0, fp+1)
SCADcoef[c(1, pickind+1)]=tempw # SCADresult$w
  
 if (detailed) {
        return(list(wt.initsoln = wt.initsoln, SCADcoef = SCADcoef, SCADresult=SCADresult))
    }
    else {
        return(list(wt.initsoln = wt.initsoln, SCADcoef = SCADcoef))
    }  
}

#########################################################

getfinalLASSOcoef <- function(x, y, pickind, folds=folds, eps0=1e-3, family=binomial(), tune.method="AIC", inittype='NoPen') {

fp=ncol(x)
fn=nrow(x)
ones=rep(1, fn)


wt.initsoln=rep(0, fp+1)

SCADresult=switch(tune.method,
CV=CVscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, folds=folds, family=family),
AIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC='AIC', eps0=eps0),
BIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC='BIC', eps0=eps0))
         

tempw=SCADresult$w
tempw[abs(tempw)<eps0]=0
SCADcoef=rep(0, fp+1)
SCADcoef[c(1, pickind+1)]=tempw # SCADresult$w
  
  
return(list(wt.initsoln=wt.initsoln, SCADcoef=SCADcoef, SCADresult=SCADresult))
}
      
###############################################
GLMvanISISscad <- function(x, y, nsis=NULL, family=binomial(), folds=folds, rank.method="obj",
  eps0=1e-3, inittype='NoPen', tune.method="AIC",  ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5) {
### ISIStypeCumulative: TRUE     selected variable is     curmulated
#                        FALSE                         NOT curmulated
if((inittype!='NoPen')&&(inittype!='L1')) stop("inittype must be either 'NoPen' or 'L1'")
  x=as.matrix(x)
  fp=ncol(x)
  fn=nrow(x)
  if(is.null(nsis)) nsis=floor(min(fp,fn/log(fn)/4))
  ones=rep(1, fn)
  tempdev0=NULL
  tempmycoef=NULL
  for (i in 1:fp) {
    tempglmfit=glm.fit(cbind(ones, x[,i]), y, family=family)
    tempdev0[i]=tempglmfit$deviance
    tempmycoef[i]=coef(tempglmfit)[2]
  }
  used.dev=switch(rank.method,  obj=tempdev0, coeff=-abs(tempmycoef))
  tempdev0.sort=sort(used.dev, method = "sh", index=TRUE)
  initRANKorder= tempdev0.sort$ix
  SISind= sort(tempdev0.sort$ix[1:nsis])
  if(!DOISIS) return(list(initRANKorder=initRANKorder, SISind=SISind,  nsis=nsis))

  allind=seq(1, fp, 1)
#  pickind=sort(tempdev0.sort$ix[1:nsis])
#  pickind=pickind[1:floor(2*length(pickind)/3)]
 pickind=sort(tempdev0.sort$ix[1:floor(nsis*2/3)])

if(inittype=='L1') {
    L1result=switch(tune.method,
            CV=CVscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, folds=folds, family=family),
            AIC=AICBICscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0),
            BIC=AICBICscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0))
  wt.initsoln=L1result$w
}
if(inittype=='NoPen') {
  wt.initsoln=coef(glm.fit(cbind(ones, x[, pickind]), y, family=family))
}
   
SCADresult=switch(tune.method,
          CV=CVscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, folds=folds, family=family),
          AIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0),
          BIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0)) 
cur.coef=SCADresult$w

nonzeroindc=which(abs(cur.coef[-1])>eps0)

ISISind=sort(pickind[nonzeroindc])
test=if(length(ISISind)<nsis) 1 else 0
normal.exit=1
  
detail.pickind=NULL
detail.ISISind=NULL
curloop=1
detail.pickind=as.list(detail.pickind)
detail.ISISind=as.list(detail.ISISind)
detail.pickind[[curloop]]=pickind
detail.ISISind[[curloop]]=ISISind
  
while(test) {
    oldISISind=ISISind
    remind=setdiff(allind, ISISind)
    tempdev0=NULL
    tempmycoef=NULL
    for (i in 1:length(remind)) {
      tempglmfit=glm.fit(cbind(ones, x[,c(remind[i], ISISind)]), y, family=family)
      tempdev0[i]=tempglmfit$deviance
      tempmycoef[i]=coef(tempglmfit)[2]
     }
    used.dev=switch(rank.method,  obj=tempdev0, coeff=-abs(tempmycoef))
    tempdev0.sort=sort(used.dev, method = "sh", index=TRUE)
    new.pickind=sort(remind[tempdev0.sort$ix[1:(nsis-length(ISISind))]])
    
    pickind=c(ISISind, new.pickind)

  if(inittype=='L1') {
     L1result=switch(tune.method,
            CV=CVscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, folds=folds, family=family),
            AIC=AICBICscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0),
            BIC=AICBICscadglm(x[, pickind], y, wt.initsoln=rep(0, length(pickind)+1), dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0))
   wt.initsoln=L1result$w
  }
  if(inittype=='NoPen') {
    wt.initsoln=coef(glm.fit(cbind(ones, x[, pickind]), y, family=family))
   }
   
  if(ISIStypeCumulative) nopenalty.subset=seq(length(ISISind)) else   nopenalty.subset=NULL
  
   SCADresult=switch(tune.method,
            CV=CVscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, folds=folds, family=family, nopenalty.subset=nopenalty.subset),
            AIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0, nopenalty.subset=nopenalty.subset),
            BIC=AICBICscadglm(x[, pickind], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0, nopenalty.subset=nopenalty.subset)) 

  cur.coef=SCADresult$w
  nonzeroindc=which(abs(cur.coef[-1])>eps0)

  ISISind=sort(pickind[nonzeroindc])
 
  curloop=curloop+1

  detail.pickind[[curloop]]=pickind
  detail.ISISind[[curloop]]=ISISind

   if((setequal(oldISISind, ISISind))||(length(ISISind)>=nsis)) test=0
   if(curloop>maxloop) {
     test=0
     normal.exit=0
    }

   }
   ISISind=sort(ISISind)
   SISind=sort(SISind)
   return(list(initRANKorder=initRANKorder, detail.pickind=detail.pickind, detail.ISISind=detail.ISISind,  SISind=SISind, ISISind=ISISind,  nsis=nsis, normal.exit=normal.exit))
}

###############################################################################################################################################


GLMvarISISscad <- function(x, y, nsis=NULL, family=binomial(), folds=folds, rank.method="obj",
    tune.method="AIC", vartype="First", eps0=1e-3, inittype='NoPen',  ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5) {
### ISIStypeCumulative: TRUE     selected variable is     curmulated
#                        FALSE                         NOT curmulated
if((inittype!='NoPen')&&(inittype!='L1')) stop("inittype must be either 'regularGLM' or 'regularGLM'")
x=as.matrix(x)
fp=ncol(x)
fn=nrow(x)
if(is.null(nsis)){
    if(vartype=='First') nsis=floor(min(fp,fn/log(fn))) else nsis=floor(min(fp,fn/log(fn)/4))
}
ones=rep(1, fn)
allind=seq(1, fp, 1)
fn.half=floor(fn/2)

x1=x[1:fn.half, ]
x1=scale(x1)
y1=y[1:fn.half]
ones1=ones[1:fn.half]

x2=x[(1+fn.half):fn.half, ]
x2=scale(x2)
y2=y[(1+fn.half):fn.half]
ones2=ones[(1+fn.half):fn.half]
###########################################
tempdev01=NULL
tempmycoef01=NULL
for (i in 1:fp) {
  tempglmfit=glm.fit(cbind(ones1, x1[,i]), y1, family=family)
  tempdev01[i]=tempglmfit$deviance
  tempmycoef01[i]=coef(tempglmfit)[2]
}
used.dev1=switch(rank.method, obj=tempdev01, coeff=-abs(tempmycoef01))
tempdev01.sort=sort(used.dev1, method = "sh", index=TRUE)
initRANKorder1= tempdev01.sort$ix
old.initRANKorder1=initRANKorder1
#####################
tempdev02=NULL
tempmycoef02=NULL
for (i in 1:fp) {
  tempglmfit=glm.fit(cbind(ones2, x2[,i]), y2, family=family)
  tempdev02[i]=tempglmfit$deviance
  tempmycoef02[i]=coef(tempglmfit)[2]
}
used.dev2=switch(rank.method, obj=tempdev02, coeff=-abs(tempmycoef02))
tempdev02.sort=sort(used.dev2, method = "sh", index=TRUE)
initRANKorder2= tempdev02.sort$ix
old.initRANKorder2=initRANKorder2
###########################################
firstn=nsis
pickindall=sort(intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn]))
if(vartype=="Second") {
while(length(pickindall)<nsis) {
  firstn=firstn+1
  pickindall=sort(intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn]))
}
}
SISind=sort(pickindall)

if(!DOISIS)     return(list(initRANKorder1=initRANKorder1, initRANKorder2=initRANKorder2, SISind=SISind, nsis=nsis))


#pickindall=pickindall[1:floor(2*length(pickindall)/3)]
firstn=nsis
pickindall=sort(intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn]))
if(vartype=="Second") {
while(length(pickindall)<nsis*2/3) {
  firstn=firstn+1
  pickindall=sort(intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn]))
}
}

if(inittype=='L1') {
    L1result=switch(tune.method,
            CV=CVscadglm(x[, pickindall], y, wt.initsoln=rep(0, length(pickindall)+1), dlambda=NULL, folds=folds, family=family),
            AIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=rep(0, length(pickindall)+1), dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0),
            BIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=rep(0, length(pickindall)+1), dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0))
    wt.initsoln=L1result$w
}
if(inittype=='NoPen') {
  wt.initsoln=coef(glm.fit(cbind(ones, x[, pickindall]), y, family=family))
}
   SCADresult=switch(tune.method,
            CV=CVscadglm(x[, pickindall], y, wt.initsoln=wt.initsoln, dlambda=NULL, folds=folds, family=family),
            AIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0),
            BIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0)) 

cur.coef=SCADresult$w
nonzeroindc=which(abs(cur.coef[-1])>eps0)
ISISind=sort(pickindall[nonzeroindc])

test=if(length(ISISind)<nsis) 1 else 0
normal.exit=1

curloop=1
  detail.pickind=NULL
  detail.ISISind=NULL
detail.pickind=as.list(detail.pickind)
detail.ISISind=as.list(detail.ISISind)
  detail.pickind[[curloop]]=pickindall
  detail.ISISind[[curloop]]=ISISind
while(test) {
  oldISISind=ISISind
  remind=setdiff(allind, ISISind)
  tempdev01=NULL
  tempmycoef01=NULL
  for (i in 1:length(remind)) {
    tempglmfit=glm.fit(cbind(ones1, x1[,c(remind[i], ISISind)]), y1, family=family)
    tempdev01[i]=tempglmfit$deviance
    tempmycoef01[i]=coef(tempglmfit)[2]
  }
used.dev1=switch(rank.method, obj=tempdev01, coeff=-abs(tempmycoef01))
tempdev01.sort=sort(used.dev1, method = "sh", index=TRUE)
initRANKorder1= tempdev01.sort$ix
#####################
  tempdev02=NULL
  for (i in 1:length(remind)) {
    tempglmfit=glm.fit(cbind(ones2, x2[,c(remind[i], ISISind)]), y2, family=family)
    tempdev02[i]=tempglmfit$deviance
    tempmycoef02[i]=coef(tempglmfit)[2]
  }
used.dev2=switch(rank.method, obj=tempdev02, coeff=-abs(tempmycoef02))
tempdev02.sort=sort(used.dev2, method = "sh", index=TRUE)
initRANKorder2= tempdev02.sort$ix
###########################################
 nleft=nsis-length(ISISind)
 firstn=nleft
 new.pickindall=sort(remind[intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn])])

if(vartype=="Second") { 
 while(length(new.pickindall)<nleft) {
   firstn=firstn+1
   new.pickindall=sort(remind[intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn])])
 }
 }
 pickindall=c(ISISind, new.pickindall)
 
if(inittype=='L1') {
    L1result=switch(tune.method,
            CV=CVscadglm(x[, pickindall], y, wt.initsoln=rep(0, length(pickindall)+1), dlambda=NULL, folds=folds, family=family),
            AIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=rep(0, length(pickindall)+1), dlambda=NULL, family=family, AICBIC="AIC", eps0=eps0),
            BIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=rep(0, length(pickindall)+1), dlambda=NULL, family=family, AICBIC="BIC", eps0=eps0))
    wt.initsoln=L1result$w
}
if(inittype=='NoPen') {
  wt.initsoln=coef(glm.fit(cbind(ones, x[, pickindall]), y, family=family))
}
    if(ISIStypeCumulative) nopenalty.subset=seq(length(ISISind)) else   nopenalty.subset=NULL

SCADresult=switch(tune.method,
          CV=CVscadglm(x[, pickindall], y, wt.initsoln=wt.initsoln, dlambda=NULL, folds=folds, family=family, nopenalty.subset=nopenalty.subset),
          AIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="AIC", nopenalty.subset=nopenalty.subset, eps0=eps0),
          BIC=AICBICscadglm(x[, pickindall], y, wt.initsoln=wt.initsoln, dlambda=NULL, family=family, AICBIC="BIC", nopenalty.subset=nopenalty.subset, eps0=eps0)) 


 
 cur.coef=SCADresult$w
 nonzeroindc=which(abs(cur.coef[-1])>eps0)
 ISISind=sort(pickindall[nonzeroindc])
 curloop=curloop+1
  detail.pickind[[curloop]]=pickindall
  detail.ISISind[[curloop]]=ISISind

 if((setequal(oldISISind, ISISind))||(length(ISISind)>=nsis)) test=0
 if(curloop>maxloop) {
    test=0
    normal.exit=0
 }
}
   ISISind=sort(ISISind)
   SISind=sort(SISind)
   return(list(initRANKorder1=old.initRANKorder1, initRANKorder2=old.initRANKorder2, detail.pickind=detail.pickind, detail.ISISind=detail.ISISind,  SISind=SISind, ISISind=ISISind,  nsis=nsis, normal.exit=normal.exit))
}












