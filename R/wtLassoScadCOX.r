getfinalSCADcoefCOX <- function (x, time, status, method = "efron", pickind, folds = NULL, 
    eps0 = 1e-3, tune.method = "AIC", inittype = "NoPen", detailed = FALSE) 
{
    p = ncol(x)
    n = nrow(x)
    if (inittype == "L1") {
        L1result = CVscadcox(x[, pickind], time, status, wt.initsoln = rep(0, 
            length(pickind)), method = method, dlambda = NULL, 
            folds = folds)
        wt.initsoln = L1result$w
    }
    if (inittype == "NoPen") {
        wt.initsoln = coef(coxph(Surv(time, status) ~ x[, pickind]))
    }
    SCADresult = switch(tune.method, CV = CVscadcox(x[, pickind], 
        time, status, wt.initsoln = wt.initsoln, folds = folds, 
        dlambda = NULL, method = method, weight = NULL, offset = NULL, 
        nopenalty.subset = NULL, eps0 = eps0), AIC = AICBICscadcox(x[, 
        pickind], time, status, wt.initsoln = wt.initsoln, method = method, 
        AICBIC = "AIC", weight = NULL, offset = NULL, nopenalty.subset = NULL, 
        eps0 = eps0), BIC = AICBICscadcox(x[, pickind], time, 
        status, wt.initsoln = wt.initsoln, method = method, AICBIC = "BIC", 
        weight = NULL, offset = NULL, nopenalty.subset = NULL, 
        eps0 = eps0))
    tempw = SCADresult$w
    tempw[abs(tempw) < eps0] = 0
    SCADcoef = rep(0, p)
    SCADcoef[pickind] = tempw
    if (detailed) {
        return(list(wt.initsoln = wt.initsoln, SCADcoef = SCADcoef, SCADresult=SCADresult))
    }
    else {
        return(list(wt.initsoln = wt.initsoln, SCADcoef = SCADcoef))
    }
}




########################################################################################
 "COXvanISISscad" <- function(x, time, status,  method = "efron", nsis=NULL, folds=NULL, rank.method="obj", 
 eps0=1e-3, inittype='NoPen', tune.method="AIC", ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5) {
    fn <- nrow(x)
    fp <- ncol(x)
    tempdev <- NULL
    tempcoef <- NULL
    if (is.null(nsis)) 
        nsis = floor(min(fp,fn/log(fn)/4))
    for (i in 1:fp) {
        cox.fit <- coxph(Surv(time, status) ~ x[, i])
        tempdev[i] <- -cox.fit$loglik[2]
        tempcoef[i] <- cox.fit$coef
    }
    used.list <- switch(rank.method, obj = tempdev, coeff = -abs(tempcoef))
    used.sort <- sort(used.list, method = "sh", index = TRUE)
    initRANKorder <- used.sort$ix
    SISind <- sort(initRANKorder[1:nsis])
    if (!DOISIS) 
        return(list(initRANKorder = initRANKorder, SISind = SISind, 
            nsis = nsis))
    allind = seq(1, fp, 1)
    pick.ind = initRANKorder[1:floor(2 * nsis/3)]
    if (inittype == "L1") {
        L1.result = switch(tune.method, CV = CVscadcox(x[, pick.ind], 
            time, status, wt.initsoln = rep(0, length(pick.ind)), 
            folds = folds, dlambda = NULL, method = method, weight = NULL, 
            offset = NULL, nopenalty.subset = NULL, eps0 = eps0), 
            AIC = AICBICscadcox(x[, pick.ind], time, status, 
                wt.initsoln = rep(0, length(pick.ind)), method = method, 
                AICBIC = "AIC", weight = NULL, offset = NULL, 
                nopenalty.subset = NULL, eps0 = eps0), BIC = AICBICscadcox(x[, 
                pick.ind], time, status, wt.initsoln = rep(0, 
                length(pick.ind)), method = method, AICBIC = "BIC", 
                weight = NULL, offset = NULL, nopenalty.subset = NULL, 
                eps0 = eps0))
        wt.initsoln = L1.result$w
    }
    if (inittype == "NoPen") {
        wt.initsoln = coef(coxph(Surv(time, status) ~ x[, pick.ind]))
    }
    SCAD.result = switch(tune.method, CV = CVscadcox(x[, pick.ind], 
        time, status, wt.initsoln = wt.initsoln, folds = folds, 
        dlambda = NULL, method = method, weight = NULL, offset = NULL, 
        nopenalty.subset = NULL, eps0 = eps0), AIC = AICBICscadcox(x[, 
        pick.ind], time, status, wt.initsoln = wt.initsoln, method = method, 
        AICBIC = "AIC", weight = NULL, offset = NULL, nopenalty.subset = NULL, 
        eps0 = eps0), BIC = AICBICscadcox(x[, pick.ind], time, 
        status, wt.initsoln = wt.initsoln, method = method, AICBIC = "BIC", 
        weight = NULL, offset = NULL, nopenalty.subset = NULL, 
        eps0 = eps0))
    current.coef = SCAD.result$w
    non.zero.ind = which(abs(current.coef) > eps0)
    ISISind = sort(pick.ind[non.zero.ind])
    test = if (length(ISISind) < nsis) 
        1
    else 0
    normal.exit = 1
    detail.pickind = NULL
    detail.ISISind = NULL
    detail.pickind = as.list(detail.pickind)
    detail.ISISind = as.list(detail.ISISind)
    curloop = 1
    detail.pickind[[curloop]] = pick.ind
    detail.ISISind[[curloop]] = ISISind
    if (length(ISISind) == 0) {
        warnings("No variable was selected by this variant of ISIS.")
        return(list(initRANKorder = initRANKorder, detail.pickind = detail.pickind, 
            detail.ISISind = detail.ISISind, SISind = SISind, 
            ISISind = ISISind, nsis = nsis, normal.exit = normal.exit))
    }
    while (test) {
        oldISISind = ISISind
        remind = setdiff(allind, ISISind)
        tempdev0 = NULL
        tempmycoef = NULL
        for (i in 1:length(remind)) {
            tempcox.fit = coxph(Surv(time, status) ~ x[, c(remind[i], 
                ISISind)])
            tempdev0[i] = -tempcox.fit$loglik[2]
            tempmycoef[i] = tempcox.fit$coef[1]
        }
        used.dev = switch(rank.method, obj = tempdev0, coeff = -abs(tempmycoef))
        tempdev0.sort = sort(used.dev, method = "sh", index = TRUE)
        new.pickind = sort(remind[tempdev0.sort$ix[1:(nsis - 
            length(ISISind))]])
        pick.ind = c(ISISind, new.pickind)
        if (inittype == "L1") {
            L1.result = switch(tune.method, CV = CVscadcox(x[, 
                pick.ind], time, status, wt.initsoln = rep(0, 
                length(pick.ind)), folds = folds, dlambda = NULL, 
                method = method, weight = NULL, offset = NULL, 
                nopenalty.subset = NULL, eps0 = eps0), AIC = AICBICscadcox(x[, 
                pick.ind], time, status, wt.initsoln = rep(0, 
                length(pick.ind)), method = method, AICBIC = "AIC", 
                weight = NULL, offset = NULL, nopenalty.subset = NULL, 
                eps0 = eps0), BIC = AICBICscadcox(x[, pick.ind], 
                time, status, wt.initsoln = rep(0, length(pick.ind)), 
                method = method, AICBIC = "BIC", weight = NULL, 
                offset = NULL, nopenalty.subset = NULL, eps0 = eps0))
            wt.initsoln = L1.result$w
        }
        if (inittype == "NoPen") {
            wt.initsoln = coef(coxph(Surv(time, status) ~ x[, 
                pick.ind]))
        }
        if (ISIStypeCumulative) 
            nopenalty.subset = seq(length(ISISind))
        else nopenalty.subset = NULL
        SCAD.result = switch(tune.method, CV = CVscadcox(x[, 
            pick.ind], time, status, wt.initsoln = wt.initsoln, 
            folds = folds, dlambda = NULL, method = method, weight = NULL, 
            offset = NULL, nopenalty.subset = NULL, eps0 = eps0), 
            AIC = AICBICscadcox(x[, pick.ind], time, status, 
                wt.initsoln = wt.initsoln, method = method, AICBIC = "AIC", 
                weight = NULL, offset = NULL, nopenalty.subset = NULL, 
                eps0 = eps0), BIC = AICBICscadcox(x[, pick.ind], 
                time, status, wt.initsoln = wt.initsoln, method = method, 
                AICBIC = "BIC", weight = NULL, offset = NULL, 
                nopenalty.subset = NULL, eps0 = eps0))
        cur.coef = SCAD.result$w
        nonzeroindc = which(abs(cur.coef) > eps0)
        ISISind = sort(pick.ind[nonzeroindc])
        curloop = curloop + 1
        detail.pickind[[curloop]] = pick.ind
        detail.ISISind[[curloop]] = ISISind
        if ((setequal(oldISISind, ISISind)) || (length(ISISind) >= 
            nsis)) 
            test = 0
        if (curloop > maxloop) {
            test = 0
            normal.exit = 0
        }
    }
    ISISind = sort(ISISind)
    SISind = sort(SISind)
    return(list(initRANKorder = initRANKorder, detail.pickind = detail.pickind, 
        detail.ISISind = detail.ISISind, SISind = SISind, ISISind = ISISind, 
        nsis = nsis, normal.exit = normal.exit))
}

 
 
COXvarISISscad <- function(x, time, status,  method = "efron", nsis=NULL, folds=NULL, rank.method="obj", 
 eps0=1e-3, inittype='NoPen', tune.method="AIC", vartype="First",  ISIStypeCumulative=FALSE, DOISIS=TRUE, maxloop=5) {
    if ((inittype != "NoPen") && (inittype != "L1")) 
        stop("wrong inittype")
    x = as.matrix(x)
    fp = ncol(x)
    fn = nrow(x)
    if (is.null(nsis)) {
        if (vartype == "First") 
            nsis = floor(min(fp,fn/log(fn)))
        else nsis = floor(min(fp,fn/log(fn)/4))
    }
    allind = seq(1, fp, 1)
    fn.half = floor(fn/2)
    x1 = x[1:fn.half, ]
    x1 = scale(x1)
    time1 = time[1:fn.half]
    status1 = status[1:fn.half]
    x2 = x[(1:fn.half) + fn.half, ]
    x2 = scale(x2)
    time2 = time[(1:fn.half) + fn.half]
    status2 = status[(1:fn.half) + fn.half]
    tempdev01 = NULL
    tempmycoef01 = NULL
    for (i in 1:fp) {
        tempcoxfit <- coxph(Surv(time1, status1) ~ x1[, i])
        tempdev01[i] <- -tempcoxfit$loglik[2]
        tempmycoef01[i] <- tempcoxfit$coef
    }
    used.dev1 = switch(rank.method, obj = tempdev01, coeff = -abs(tempmycoef01))
    tempdev01.sort = sort(used.dev1, method = "sh", index = TRUE)
    initRANKorder1 = tempdev01.sort$ix
    old.initRANKorder1 = initRANKorder1
    tempdev02 = NULL
    tempmycoef02 = NULL
    for (i in 1:fp) {
        tempcoxfit <- coxph(Surv(time2, status2) ~ x2[, i])
        tempdev02[i] <- -tempcoxfit$loglik[2]
        tempmycoef02[i] <- tempcoxfit$coef
    }
    used.dev2 = switch(rank.method, obj = tempdev02, coeff = -abs(tempmycoef02))
    tempdev02.sort = sort(used.dev2, method = "sh", index = TRUE)
    initRANKorder2 = tempdev02.sort$ix
    old.initRANKorder2 = initRANKorder2
    firstn = nsis
    pickindall = intersect(tempdev01.sort$ix[1:firstn], tempdev02.sort$ix[1:firstn])
    if (vartype == "Second") {
        while (length(pickindall) < nsis) {
            firstn = firstn + 1
            pickindall = intersect(tempdev01.sort$ix[1:firstn], 
                tempdev02.sort$ix[1:firstn])
        }
    }
    SISind = sort(pickindall)
    if (length(pickindall) == 0) {
        warning("no intersection between the two samples")
        return(list(initRANKorder1 = initRANKorder1, initRANKorder2 = initRANKorder2, 
            SISind = SISind, nsis = nsis, normal.exit = 0))
    }
    if (!DOISIS) 
        return(list(initRANKorder1 = initRANKorder1, initRANKorder2 = initRANKorder2, 
            SISind = SISind, nsis = nsis))
    if (length(pickindall) > 2) 
        pickindall = pickindall[1:floor(2 * length(pickindall)/3)]
    if (length(pickindall) == 0) {
        warning("no intersection between the two samples")
    }
    if (inittype == "L1") {
        L1.result = switch(tune.method, CV = CVscadcox(x[, pickindall], 
            time, status, wt.initsoln = rep(0, length(pickindall)), 
            folds = folds, dlambda = NULL, method = method, weight = NULL, 
            offset = NULL, nopenalty.subset = NULL, eps0 = eps0), 
            AIC = AICBICscadcox(x[, pickindall], time, status, 
                wt.initsoln = rep(0, length(pickindall)), method = method, 
                AICBIC = "AIC", weight = NULL, offset = NULL, 
                nopenalty.subset = NULL, eps0 = eps0), BIC = AICBICscadcox(x[, 
                pickindall], time, status, wt.initsoln = rep(0, 
                length(pickindall)), method = method, AICBIC = "BIC", 
                weight = NULL, offset = NULL, nopenalty.subset = NULL, 
                eps0 = eps0))
        wt.initsoln = L1.result$w
    }
    if (inittype == "NoPen") {
        wt.initsoln = coef(coxph(Surv(time, status) ~ x[, pickindall]))
    }
    SCAD.result = switch(tune.method, CV = CVscadcox(x[, pickindall], 
        time, status, wt.initsoln = wt.initsoln, folds = folds, 
        dlambda = NULL, method = method, weight = NULL, offset = NULL, 
        nopenalty.subset = NULL, eps0 = eps0), AIC = AICBICscadcox(x[, 
        pickindall], time, status, wt.initsoln = wt.initsoln, 
        method = method, AICBIC = "AIC", weight = NULL, offset = NULL, 
        nopenalty.subset = NULL, eps0 = eps0), BIC = AICBICscadcox(x[, 
        pickindall], time, status, wt.initsoln = wt.initsoln, 
        method = method, AICBIC = "BIC", weight = NULL, offset = NULL, 
        nopenalty.subset = NULL, eps0 = eps0))
    cur.coef = SCAD.result$w
    nonzeroindc = which(abs(cur.coef) > eps0)
    ISISind = sort(pickindall[nonzeroindc])
    test = if (length(ISISind) < nsis) 
        1
    else 0
    normal.exit = 1
    curloop = 1
    detail.pickind = NULL
    detail.ISISind = NULL
    detail.pickind = as.list(detail.pickind)
    detail.ISISind = as.list(detail.ISISind)
    detail.pickind[[curloop]] = pickindall
    detail.ISISind[[curloop]] = ISISind
    if (length(ISISind) == 0) {
        warnings("No variable was selected by this variant of ISIS.")
        normal.exit = 0
        return(list(initRANKorder1 = old.initRANKorder1, initRANKorder2 = old.initRANKorder2, 
            detail.pickind = detail.pickind, detail.ISISind = detail.ISISind, 
            SISind = SISind, ISISind = ISISind, nsis = nsis, 
            normal.exit = normal.exit))
    }
    while (test) {
        oldISISind = ISISind
        remind = setdiff(allind, ISISind)
        tempdev01 = NULL
        tempmycoef01 = NULL
        for (i in 1:length(remind)) {
            tempcoxfit <- coxph(Surv(time1, status1) ~ x1[, c(remind[i], 
                ISISind)])
            tempdev01[i] <- -tempcoxfit$loglik[2]
            tempmycoef01[i] <- tempcoxfit$coef[1]
        }
        used.dev1 = switch(rank.method, obj = tempdev01, coeff = -abs(tempmycoef01))
        tempdev01.sort = sort(used.dev1, method = "sh", index = TRUE)
        initRANKorder1 = tempdev01.sort$ix
        tempdev02 = NULL
        tempmycoef02 = NULL
        for (i in 1:length(remind)) {
            tempcoxfit <- coxph(Surv(time2, status2) ~ x2[, c(remind[i], 
                ISISind)])
            tempdev02[i] <- -tempcoxfit$loglik[2]
            tempmycoef02[i] <- tempcoxfit$coef[1]
        }
        used.dev2 = switch(rank.method, obj = tempdev02, coeff = -abs(tempmycoef02))
        tempdev02.sort = sort(used.dev2, method = "sh", index = TRUE)
        initRANKorder2 = tempdev02.sort$ix
        nleft <- nsis - length(oldISISind)
        firstn <- nleft
        new.pickindall = sort(remind[intersect(tempdev01.sort$ix[1:firstn], 
            tempdev02.sort$ix[1:firstn])])
        if (vartype == "Second") {
            while (length(new.pickindall) < nleft) {
                firstn = firstn + 1
                new.pickindall = sort(remind[intersect(tempdev01.sort$ix[1:firstn], 
                  tempdev02.sort$ix[1:firstn])])
            }
        }
        pickindall = c(ISISind, new.pickindall)
        if (length(pickindall) == 0) {
            warning("no intersection between the two samples")
        }
        if (inittype == "L1") {
            L1.result = switch(tune.method, CV = CVscadcox(x[, 
                pickindall], time, status, wt.initsoln = rep(0, 
                length(pickindall)), folds = folds, dlambda = NULL, 
                method = method, weight = NULL, offset = NULL, 
                nopenalty.subset = NULL, eps0 = eps0), AIC = AICBICscadcox(x[, 
                pickindall], time, status, wt.initsoln = rep(0, 
                length(pickindall)), method = method, AICBIC = "AIC", 
                weight = NULL, offset = NULL, nopenalty.subset = NULL, 
                eps0 = eps0), BIC = AICBICscadcox(x[, pickindall], 
                time, status, wt.initsoln = rep(0, length(pickindall)), 
                method = method, AICBIC = "BIC", weight = NULL, 
                offset = NULL, nopenalty.subset = NULL, eps0 = eps0))
            wt.initsoln = L1.result$w
        }
        if (inittype == "NoPen") {
            wt.initsoln = coef(coxph(Surv(time, status) ~ x[, 
                pickindall]))
        }
        if (ISIStypeCumulative) 
            nopenalty.subset = seq(length(ISISind))
        else nopenalty.subset = NULL
        SCAD.result = switch(tune.method, CV = CVscadcox(x[, 
            pickindall], time, status, wt.initsoln = wt.initsoln, 
            folds = folds, dlambda = NULL, method = method, weight = NULL, 
            offset = NULL, nopenalty.subset = nopenalty.subset, 
            eps0 = eps0), AIC = AICBICscadcox(x[, pickindall], 
            time, status, wt.initsoln = wt.initsoln, method = method, 
            AICBIC = "AIC", weight = NULL, offset = NULL, nopenalty.subset = nopenalty.subset, 
            eps0 = eps0), BIC = AICBICscadcox(x[, pickindall], 
            time, status, wt.initsoln = wt.initsoln, method = method, 
            AICBIC = "BIC", weight = NULL, offset = NULL, nopenalty.subset = nopenalty.subset, 
            eps0 = eps0))
        cur.coef = SCAD.result$w
        nonzeroindc = which(abs(cur.coef) > eps0)
        ISISind = sort(pickindall[nonzeroindc])
        curloop = curloop + 1
        detail.pickind[[curloop]] = pickindall
        detail.ISISind[[curloop]] = ISISind
        if ((setequal(oldISISind, ISISind)) || (length(ISISind) >= 
            nsis)) 
            test = 0
        if (curloop > maxloop) {
            warnings("maximum loop achieved...")
            test = 0
            normal.exit = 0
        }
    }
    ISISind = sort(ISISind)
    SISind = sort(SISind)
    return(list(initRANKorder1 = old.initRANKorder1, initRANKorder2 = old.initRANKorder2, 
        detail.pickind = detail.pickind, detail.ISISind = detail.ISISind, 
        SISind = SISind, ISISind = ISISind, nsis = nsis, normal.exit = normal.exit))
}

 
 

############################################################################## 
"CVscadcox" <- function(x, time, status, 
wt.initsoln=NULL, folds=NULL, dlambda=NULL, 
method = "efron", lambda2=0,weight=NULL, offset=NULL, 
nopenalty.subset=NULL,eps0=1e-16) { 
x=as.matrix(x) 
p=ncol(x)
n=nrow(x) 
if(is.null(offset)) offset=rep(0, nrow (x)) 
if(is.null(weight)) weight=rep(1, nrow (x)) 
if(is.null(wt.initsoln)) wt.initsoln=rep(0, p)


if(is.null(folds)) {
temp= sample(1:n, n, replace = FALSE)
kfold=10
for(i in 1:kfold){
  folds[[i]]=setdiff(1:n, temp[seq(i, n, kfold)])
  }
}

max.lambda = findSCADMaxLambda.cox(x, time, status, wt.initsoln) 

if(max.lambda<4){
  if(is.null(dlambda))dlambda=max.lambda/40
  lambda.cand=seq(max.lambda,0, -dlambda)
  lambda.cand=lambda.cand[lambda.cand>0]
  }
if(max.lambda>=4){
  lambda.cand=max.lambda*(2^(seq(0, -10, -0.1)))
  }
  
n.folds=length(folds)
nobs.fold=NULL
for (i in 1:n.folds) nobs.fold[i]=length(folds[[i]])
  
old.soln=NULL
old.soln.save=NULL
tuneer=NULL
iloop=1
tuneer[iloop]=0

    for (j in 1:n.folds){
    tmptime=time[folds[[j]]]
    tmpstatus=status[folds[[j]]]
    tmptimeremind=time[-folds[[j]]]
    tmpstatusremind=status[-folds[[j]]]
    }


    
for (j in 1:n.folds) {
  old.soln[[j]]=scadcox(x[folds[[j]],], time[folds[[j]]], status[folds[[j]]], method, wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=NULL, weight = weight[folds[[j]]], lambda2, nopenalty.subset=nopenalty.subset)
  tuneer[iloop]=tuneer[[iloop]]-logplik(x[-folds[[j]],], time[-folds[[j]]], status[-folds[[j]]],  old.soln[[j]], method)
}
old.soln.save[[iloop]]=old.soln[[1]]

temp.initsoln=wt.initsoln

for (iloop in 2:length(lambda.cand)) {
  if((max(abs(temp.initsoln))>0) && (min(abs(temp.initsoln[abs(temp.initsoln)>0]))>3.7*lambda.cand[iloop-1])) break
  tuneer[iloop]=0
  for (j in 1:n.folds) {
    old.soln[[j]]=scadcox(x[folds[[j]],], time[folds[[j]]], status[folds[[j]]], method, wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=old.soln[[j]], weight = weight[folds[[j]]], lambda2, nopenalty.subset=nopenalty.subset)
    tuneer[iloop]=tuneer[[iloop]]-logplik(x[-folds[[j]],], time[-folds[[j]]], status[-folds[[j]]],  old.soln[[j]], method)
  }
old.soln.save[[iloop]]=old.soln[[1]]

 }
best.lambda.ind=which.min(tuneer) 
w=scadcox(x, time, status, method, wt.initsoln=wt.initsoln, lambda=lambda.cand[best.lambda.ind], initsoln=old.soln.save[[best.lambda.ind]], weight = weight, nopenalty.subset=nopenalty.subset)
return(list(wt.initsoln=wt.initsoln, tuneer=tuneer, best.lambda.ind=best.lambda.ind, max.lambda=max.lambda, w=w, used.lambda=lambda.cand[1:iloop], best.lambda=lambda.cand[best.lambda.ind]))
} 


###################################################################################################
############################################################################## 
"CVlassocox" <- function(x, time, status, 
wt.initsoln=NULL, folds=NULL, dlambda=NULL, 
method = "efron", lambda2=0,weight=NULL, offset=NULL, 
nopenalty.subset=NULL,eps0=1e-16) { 
x=as.matrix(x) 
n=nrow(x) 
p=ncol(x)
if(is.null(offset)) offset=rep(0, nrow (x)) 
if(is.null(weight)) weight=rep(1, nrow (x)) 
if(is.null(wt.initsoln)) wt.initsoln=rep(0, p)


if(is.null(folds)) {
temp= sample(1:n, n, replace = FALSE)
kfold=10
for(i in 1:kfold){
  folds[[i]]=setdiff(1:n, temp[seq(i, n, kfold)])
  }
}

max.lambda = findLASSOMaxLambda.cox(x, time, status) 

if(max.lambda<4){
  if(is.null(dlambda))dlambda=max.lambda/40
  lambda.cand=seq(max.lambda,0, -dlambda)
  lambda.cand=lambda.cand[lambda.cand>0]
  }
if(max.lambda>=4){
  lambda.cand=max.lambda*(2^(seq(0, -10, -0.1)))
  }
  
n.folds=length(folds)
nobs.fold=NULL
for (i in 1:n.folds) nobs.fold[i]=length(folds[[i]])
  
old.soln=NULL
old.soln.save=NULL
tuneer=NULL
iloop=1
tuneer[iloop]=0

#    for (j in 1:n.folds){
#    tmptime=time[folds[[j]]]
#    tmpstatus=status[folds[[j]]]
#    tmptimeremind=time[-folds[[j]]]
#    tmpstatusremind=status[-folds[[j]]]
#    }


    
for (j in 1:n.folds) {
  old.soln[[j]]=wtlassocox(x[folds[[j]],], time[folds[[j]]], status[folds[[j]]], method, lassoweight=rep(n*lambda.cand[iloop],p), initsoln=NULL, weight = weight[folds[[j]]], lambda2)$w
  tuneer[iloop]=tuneer[iloop]-logplik(x[-folds[[j]],], time[-folds[[j]]], status[-folds[[j]]],  old.soln[[j]], method)
}
old.soln.save[[iloop]]=old.soln[[1]]

temp.initsoln=wt.initsoln

for (iloop in 2:length(lambda.cand)) {
  print(iloop)
  print(sum(abs(old.soln.save[[iloop-1]])>1e-3))
  #if((max(abs(temp.initsoln))>0) && (min(abs(temp.initsoln[abs(temp.initsoln)>0]))>3.7*lambda.cand[iloop-1])) break
  tuneer[iloop]=0
  for (j in 1:n.folds) {
    old.soln[[j]]=wtlassocox(x[folds[[j]],], time[folds[[j]]], status[folds[[j]]], method, lassoweight=rep(n*lambda.cand[iloop],p), initsoln=old.soln[[j]], weight = weight[folds[[j]]], lambda2)$w
    tuneer[iloop]=tuneer[iloop]-logplik(x[-folds[[j]],], time[-folds[[j]]], status[-folds[[j]]],  old.soln[[j]], method)
  }
old.soln.save[[iloop]]=old.soln[[1]]

 }
best.lambda.ind=which.min(tuneer) # note: this selection of tuning parameter may not be very good, an alternative one is given below
w=wtlassocox(x, time, status, method, lassoweight=rep(n*lambda.cand[best.lambda.ind],p), initsoln=old.soln.save[[best.lambda.ind]], weight = weight)$w
#return(list(x=x, y=y, wt.initsoln=wt.initsoln, tuneer=tuneer, best.lambda.ind=best.lambda.ind, maxlambda=maxlambda, w=w, used.lambda=lambda.cand[1:iloop]))
return(list(wt.initsoln=wt.initsoln, tuneer=tuneer, best.lambda.ind=best.lambda.ind, max.lambda=max.lambda, w=w, used.lambda=lambda.cand[1:iloop], best.lambda=lambda.cand[best.lambda.ind]))
} 

############################################################################################################################################################ 
"AICBICscadcox" <- function(x, time, status, 
wt.initsoln=NULL, dlambda=NULL,  method = 
"efron",lambda2=0, AICBIC='AIC', weight=NULL, offset=NULL, 
nopenalty.subset=NULL, eps0=1e-16) { 

x=as.matrix(x) 
if(is.null(offset)) offset=rep(0, nrow(x)) 
if(is.null(weight)) weight=rep(1, nrow(x)) 
if(is.null(wt.initsoln)) wt.initsoln=rep(0, ncol(x))

max.lambda = findSCADMaxLambda.cox(x, time, status,  wt.initsoln) 

if(max.lambda<4) {
  if(is.null(dlambda)) dlambda=max.lambda/40
  lambda.cand=seq(max.lambda, 0, -dlambda)
  lambda.cand=lambda.cand[lambda.cand>0]
 } 
 if(max.lambda>=4)  {
  lambda.cand=max.lambda*(2^(seq(0, -10, -0.1)))
 }
 
AB.factor=switch(AICBIC, AIC=2, BIC=log(nrow(x)))

old.soln=NULL

tuneer=NULL
iloop=1

old.soln[[iloop]]=scadcox(x, time, status, method, wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=NULL, weight = weight, lambda2, nopenalty.subset=nopenalty.subset)
tuneer[iloop]=length(which(abs(old.soln[[iloop]])>eps0))*AB.factor-logplik(x, time, status,   old.soln[[iloop]], method)


for (iloop in 2:length(lambda.cand)) {
old.soln[[iloop]]=scadcox(x, time, status,  method, wt.initsoln=wt.initsoln, lambda=lambda.cand[iloop], initsoln=old.soln[[iloop-1]], weight = weight, lambda2, nopenalty.subset=nopenalty.subset)
tuneer[iloop]=length(which(abs(old.soln[[iloop]])>eps0))*AB.factor-logplik(x, time, status, old.soln[[iloop]], method)
}

best.lambda.ind=which.min(tuneer) # note: this selection of tuning parameter may not be very good, an alternative one is given below

w=scadcox(x, time, status, method, wt.initsoln=wt.initsoln, lambda=lambda.cand[best.lambda.ind], initsoln=old.soln[[best.lambda.ind]], weight = weight, nopenalty.subset=nopenalty.subset)
return(list(wt.initsoln=wt.initsoln, tuneer=tuneer, best.lambda.ind=best.lambda.ind, max.lambda=max.lambda, w=w, used.lambda=lambda.cand[1:iloop], best.lambda=lambda.cand[best.lambda.ind]))
}    # end of CVscadcox function 
  
##############################################################################
"scadcox" <- function(x, time, status, method="efron", wt.initsoln=NULL, lambda, 
    initsoln=NULL, weight = NULL, function.precision=1e-8, nopenalty.subset=NULL){
    x = as.matrix(x)
    if (is.null(weight)) 
        weight = rep(1, nrow(x))
    if (is.null(wt.initsoln)) 
        wt.initsoln = rep(0, ncol(x))
    n <- nrow(x)
    p <- ncol(x)
    lasso.weight <- n * scadderiv(abs(wt.initsoln), 3.7, lambda)
    lasso.weight[nopenalty.subset] <- 0
    return(wtlassocox(x, time, status, method, lasso.weight, 
        initsoln, weight, lambda2 = 0, function.precision)$w)
}


##############################################################################
"fullscadcox" <- function(x, time, status, method="efron", lambda, 
    initsoln=NULL, weight = NULL, function.precision=1e-8, nopenalty.subset=NULL, eps0=1e-6){

 x = as.matrix(x)
    if (is.null(weight)) 
        weight = rep(1, nrow(x))
    if (is.null(initsoln)) 
        initsoln = rep(0, ncol(x))
    n <- nrow(x)
    p <- ncol(x)
    wold = initsoln
    test = 1
    while (test) {
        lasso.weight <- n * scadderiv(abs(wold), 3.7, lambda)
        lasso.weight[nopenalty.subset] <- 0
        wnew = wtlassocox(x, time, status, method, lasso.weight, 
            initsoln, weight, lambda2 = 0, function.precision)$w
        if (sum((wold - wnew)^2) < eps0) 
            test = 0
        wold = wnew
    }
    return(wnew)
}

###############################################################################
"findSCADMaxLambda.cox" <- function(x, y, d, wt.initsoln=wt.initsoln) {
  x=as.matrix(x)
    eta <- rep(0,ncol(x)) 
    n <- nrow(x)
    for (i in 1:length(d)){
        if(d[i]==1){
        riskset=(y>=y[i])
        R=sum(riskset)
        if(ncol(x)==1) xrisksetsum=sum(x[riskset,])    else if(R==1) 
         xrisksetsum=x[riskset,]      else
         xrisksetsum=apply(x[riskset,],2,sum)
        eta <- eta-(x[i,]-1/R*xrisksetsum)
        }
        }
    
lambda=max(apply(cbind(abs(eta)/n, abs(wt.initsoln)), 1, scadafindmaxlam))
return(lambda)
 }

###############################################################################
"findLASSOMaxLambda.cox" <- function(x, y, d) {
  x=as.matrix(x)
    eta <- rep(0,ncol(x)) 
    n <- nrow(x)
    for (i in 1:length(d)){
        if(d[i]==1){
        riskset=(y>=y[i])
        R=sum(riskset)
        if(ncol(x)==1) xrisksetsum=sum(x[riskset,])    else if(R==1) 
         xrisksetsum=x[riskset,]      else
         xrisksetsum=apply(x[riskset,],2,sum)
        eta <- eta-(x[i,]-1/R*xrisksetsum)
        }
        }
    
lambda=max(cbind(abs(eta)/n))
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
##############################################################################
scadderiv <- function(ftheta, fa, flambda) {
 return(flambda*(1-(1-apply(as.matrix(fa*flambda-abs(ftheta)), 1, positivepart)/((fa-1)*flambda))*as.numeric(abs(ftheta)>flambda)))
}
############################################################################################################################################################
positivepart <- function(fx) {
  return(ifelse(fx>=0, fx, 0))
}
############################################################################################################################################################
"logplik" <- function(x, time, status, b, method = c("breslow", "efron"), return.all=FALSE)
  {
    x=as.matrix(x)
    method <- match.arg(method)
    n <- length(time)
    o <- order(status, decreasing=T)
    oo <- o[order(time[o])]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n]==time[i] & status[i:n]==1)
    complete <- which(status==1)
    nnc <- length(complete)
    if (nnc == 0) {
      stop("No complete observation found. Failed to compute partial likelihood.")
    }
    dmat <- matrix(0, n, nnc)
    for (i in 1:nnc) {
      dmat[time >= time[complete[i]], i] <- 1
      if (method=="efron") {
        if (rept[complete[i]] > 0) {
          tie <- time==time[complete[i]] & status==1
          di <- max(rept[tie])
          dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]])/di
        }
      }
    }
    eta <- x %*% b
    eeta <- exp(eta)
    k <- ncol(eta)
    loglik <- rep(0, k)
    for (i in 1:k) {
      w <- dmat * eeta[oo, i]
      wsum <- apply(w, 2, sum)
      loglik[i] <- sum(eta[oo, i][status==1]) - sum(log(wsum))
    }
    if (return.all)
      return(list(loglik=loglik, w=scale(w, F, wsum), eta=eta, dmat=dmat, oo=oo))
    else return(loglik)    
  }



##############################################################################
getRept <- function(time,status){
    n <- length(time)
    o <- order(status)
    oo <- o[order(time[o], decreasing=T)]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n]==time[i] & status[i:n]==1)
    return(list(rept=rept, time=time, status=status, oo=oo))
    }



##############################################################################
"wtlassocox" <- function(x, time, status, method="efron", lassoweight=NULL, initsoln=NULL,  weight = NULL, lambda2=0, function.precision=1e-8) {

    call <- match.call()
    x = as.matrix(x)
    nobs <- nrow(x)
    m <- ncol(x)
    o <- order(status)
    oo <- o[order(time[o], decreasing = T)]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, nobs)
    for (i in 1:nobs) rept[i] <- sum(time[i:nobs] == time[i] & 
        status[i:nobs] == 1)
    x <- x[oo, ]
    x <- as.matrix(x)
    mthd = switch(method, breslow = 1, efron = 2)
    if (is.null(lassoweight)) 
        lassoweight = rep(1/nobs, m)
    if (is.null(initsoln)) 
        initsoln = rep(0, m)
    if (is.null(weight)) 
        weight = rep(1, nrow(x))
    if (length(weight) != nobs) {
        stop("Length of the weight vector != sample size")
    }
    if (any(weight < 0)) {
        stop("Negative weights are not allowed.")
    }
    if (is.null(dimnames(x)[[2]])) {
        xnames <- paste("x", seq(m), sep = "")
    }
    else xnames <- dimnames(x)[[2]]
    lambda = 0
    b0 = initsoln
    tmpa = 1:ncol(x)
    force.active = NULL
    if (is.null(force.active)) 
        a0 = b0
    else a0 = b0[-force.active]
    p <- length(tmpa)
    fk <- length(force.active)
    k <- p - fk
    if (fk > 0) 
        lassoweight <- lassoweight[-force.active]
    if (fk == 0) 
        param <- c(b0[tmpa], a0, 0, -b0[tmpa] - a0, b0[tmpa] - 
            a0)
    else param <- c(b0[tmpa], a0, 0, -b0[tmpa[-force.active]] - 
        a0, b0[tmpa[-force.active]] - a0)
    xa <- x[, tmpa, drop = FALSE]
    xstate <- rep(2, p + 3 * k + 1)
    xstate[param == 0] <- 0
    lenz <- 10 + (p + 5) * nobs + k
    zsmall <- rep(0, lenz)
    zsmall[1:7] <- c(nobs, lambda, lambda2/nobs, mthd, k, 0, 
        function.precision)
    zsmall[10 + seq((p + 3) * nobs)] <- c(as.vector(xa), time, 
        status, rept)
    zsmall[10 + (p + 5) * nobs + seq(k)] = lassoweight/nobs
    ptm <- proc.time()
    sol <- .Fortran("coxsolution", k = as.integer(k), n = as.integer(p + 
        k), nb = as.integer(p + 3 * k + 1), ne = as.integer(p + 
        5 * k), hs = as.integer(xstate), xn = as.double(param), 
        zsmall = as.double(zsmall), lenz = as.integer(lenz), 
        inform = as.integer(0))
    proc.time() - ptm
    b0[tmpa] <- sol$xn[1:p]
    if (k > 0) 
        a0 <- sol$xn[(p + 1):(p + k)]
    if (sol$inform != 0) 
        warning("\nconvergence warning in wtlassocox.\n")
    object <- list(call = call, lambda2 = lambda2, xnames = xnames, 
        weight = weight, lassoweight = lassoweight, initsoln = initsoln, 
        w = b0)
    object
}


