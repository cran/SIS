SIS <- function(x=NULL, y=NULL, family=NULL, time=NULL, status=NULL, method='efron', vartype=0, nsis=NULL, rank.method='obj', eps0=1e-3,
inittype='NoPen', tune.method='BIC', folds=NULL, post.tune.method='CV',post.tune.folds=NULL, DOISIS=TRUE, ISIStypeCumulative=FALSE, maxloop=5,
 xtune=NULL, ytune=NULL, detail=FALSE){
    if (is.null(x)) 
        stop("The covariates are missing!")
    if (!is.null(family)) {
        if (is.null(y)) 
            stop("The response is missing in Generalized Linear Model!")
        if (family$family == "binomail" & (is.null(xtune) | is.null(ytune))) 
            stop("Independent tuning data required for logit link.")
    }
    else {
        if (is.null(time) | is.null(status)) 
            stop("The data is not complete in cox model.")
    }
    n = nrow(x)
    if (is.null(nsis)) {
        if (vartype == 1) 
            nsis = floor(n/log(n))
        else nsis = floor(n/4/log(n))
    }
    if (is.null(post.tune.folds) & post.tune.method == "CV") {
        temp = sample(1:n, n, replace = FALSE)
        kfold = 10
        post.tune.folds = NULL
        for (i in 1:kfold) {
            post.tune.folds[[i]] = setdiff(1:n, temp[seq(i, n, 
                kfold)])
        }
    }
    if (!is.null(family)) {
        if (vartype == 0) {
            SISresult = GLMvanISISscad(x = x, y = y, nsis = nsis, 
                family = family, folds = folds, rank.method = rank.method, 
                eps0 = eps0, inittype = inittype, tune.method = tune.method, 
                ISIStypeCumulative = ISIStypeCumulative, DOISIS = DOISIS, 
                maxloop = maxloop)
        }
        else {
        if (vartype == 1) 
                varchar = "First"
            else varchar = "Second"
            SISresult = GLMvarISISscad(x = x, y = y, nsis = nsis, 
                family = family, folds = folds, vartype = varchar, 
                rank.method = rank.method, eps0 = eps0, inittype = inittype, 
                tune.method = tune.method, ISIStypeCumulative = ISIStypeCumulative, 
                DOISIS = DOISIS, maxloop = maxloop)
        }
        if (family$family == "binomial") {
            SIScoef = (INDEPgetfinalSCADcoef(x = x, y = y, 
                pickind = SISresult$SISind, xtune = xtune, ytune = ytune, 
                family = family, inittype = inittype))
           if (DOISIS) {
            if (length(SISresult$ISIS) == 0) {
                ISIScoef = NULL
            } else {
            ISIScoef = (INDEPgetfinalSCADcoef(x = x, y = y, 
                pickind = SISresult$ISISind, xtune = xtune, ytune = ytune, 
                family = family, inittype = inittype))
                }
        }
        }
        else {
            SIScoef = (getfinalSCADcoef(x = x, y = y, pickind = SISresult$SISind, 
                folds = post.tune.folds, eps0 = eps0, family = family, 
                tune.method = post.tune.method, inittype = "NoPen"))
                if (DOISIS) {
            if (length(SISresult$ISIS) == 0) {
                ISIScoef = NULL
            } else {
            ISIScoef = (getfinalSCADcoef(x = x, y = y, pickind = SISresult$ISISind, 
                folds = post.tune.folds, eps0 = eps0, family = family, 
                tune.method = post.tune.method, inittype = "NoPen"))
                }
        }
    }
    }
    else {
        if (vartype == 0) {
            SISresult <- COXvanISISscad(x = x, time = time, method = method, 
                folds = folds, status = status, nsis = nsis, 
                rank.method = rank.method, eps0 = eps0, inittype = inittype, 
                tune.method = tune.method, DOISIS = DOISIS, maxloop = maxloop)
        }
        else {
            if (vartype == 1) 
                varchar = "First"
            else varchar = "Second"
            SISresult <- COXvarISISscad(x = x, time = time, method = method, 
                folds = folds, status = status, nsis = nsis, 
                rank.method = rank.method, eps0 = eps0, inittype = inittype, 
                vartype = varchar, ISIStypeCumulative = ISIStypeCumulative, 
                tune.method = tune.method, DOISIS = DOISIS, maxloop = maxloop)
        }
        SIScoef = getfinalSCADcoefCOX(x = x, time = time, status = status, 
            pickind = SISresult$SIS, folds = post.tune.folds, 
            eps0 = eps0, tune.method = post.tune.method, inittype = inittype, 
            method = method)
        if (DOISIS) {
            if (length(SISresult$ISIS) == 0) {
                ISIScoef = NULL
            }
            else {
                ISIScoef = getfinalSCADcoefCOX(x = x, time = time, 
                  status = status, pickind = SISresult$ISIS, 
                  folds = post.tune.folds, eps0 = eps0, tune.method = post.tune.method, 
                  inittype = inittype, method = method)
            }
        }
    }
    if (detail == FALSE) {
        return(list(SISind = SISresult$SIS, ISISind = SISresult$ISIS, 
            SIScoef = SIScoef$SCADcoef, ISIScoef = ISIScoef$SCADcoef))
    }
    else {
        return(list(SISresult = SISresult, SIScoef = SIScoef, 
            ISIScoef = ISIScoef))
    }
}

