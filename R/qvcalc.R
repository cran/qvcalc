##  The qvcalc generic and default method.  And methods for print, summary and plot,
##  for a resulting "qv" object.


qvcalc <- function (object, ...) UseMethod("qvcalc")

qvcalc.default <- function(object, factorname = NULL, coef.indices = NULL,
                   labels = NULL, dispersion = NULL,
                   estimates = NULL, modelcall = NULL, ...)
{
    covmat <- object
    if (!is.null(labels)) rownames(covmat) <- colnames(covmat) <- labels
    n <- dim(covmat)[1]
    if (n <= 2) stop(
                    "qvcalc works only for factors with 3 or more levels")
    simple.contrasts <- function(n, levelnames = 1:n){
        result <- list()
        for (i in 1:(n-1)){
            for (j in (i+1):n){
                result[[paste(levelnames[i],
                              levelnames[j],
                              sep = ",")]] <- c(i, j)}}
        result
    }
    qvdesign <- function(n){
        nrows <- choose(n, 2)
        m <- matrix(0, nrows, n)
        indices <- simple.contrasts(n)
        for (i in 1:nrows){
            m[i, indices[[i]][1]] <- 1
            m[i, indices[[i]][2]] <- 1}
        m
    }
    level <- qvdesign(n)
    contrast.variance <- function(contrast, covmat){
        if (!(is.matrix(covmat) &&
              (dim(covmat)[1] == dim(covmat)[2])))
            stop("covmat must be a square matrix")
        n <- dim(covmat)[1]
        if (length(contrast) == n && sum(contrast) == 0)
            ## arbitrary contrast vector
            return(as.vector(contrast %*% covmat %*% contrast))
        if (length(contrast) == 2 && all(contrast %in% 1:n)){
            ## simple contrast specified as an index pair
            i <- contrast[1]
            j <- contrast[2]
            return(covmat[i,i] + covmat[j,j] - 2*covmat[i,j])}
        else stop("invalid contrast")
    }
    simple.contrast.variances <- function(n, covmat){
        if (!is.null(rownames(covmat)))
            levelnames <- rownames(covmat)
        else levelnames <- 1:n
        sapply(simple.contrasts(n, levelnames),
               function(contrast){contrast.variance(contrast, covmat)})
    }
    response <- simple.contrast.variances(n, covmat)
    if (any(response <= 0)) {
        stop("not all contrasts have positive variance")
    } else response <- log(response)
    expLinear <- structure(list(
        family = "expLinear",
        link = "exp",
        linkfun = function(mu) exp(mu),
        linkinv = function(eta) log(eta),
        variance = function(mu) rep(1, length(mu)),
        dev.resids = function(y, mu, wt) wt *
                                         ((y - mu)^2),
        aic = function(y, n, mu, wt, dev) sum(wt) *
                                          (log(dev/sum(wt) * 2 * pi) + 1) + 2,
        mu.eta = function (eta) 1/eta,
        initialize = expression({
            n <- rep(1, nobs)
            mustart <- y}),
        validmu = function(mu) TRUE),
        class = "family")
    model <- glm(response ~ 0 + level, family = expLinear)
    qv <- coef(model)
    NAs <- rep(NA, length(qv))
    if (!is.null(rownames(covmat))) names(qv) <- rownames(covmat)
    frame <- data.frame(estimate = NAs,
                        SE = sqrt(diag(covmat)),
                        quasiSE = sqrt(qv),
                        quasiVar = qv,
                        row.names = names(qv))
    if (!is.null(estimates)) frame$estimate <- estimates
    relerrs <-  sqrt(exp(- residuals(model))) - 1
    ##  The above formula was corrected in v0.8-9; it
    ##  previously said 1 - sqrt(exp(residuals(model)), which is
    ##  not what should be expected for "relative error" here.
    ##  This corrected version agrees with the Biometrika paper.
    ##  Thanks to Shaun Killingbeck for spotting this error in the
    ##  previous version.
    names(relerrs) <- names(response)
    return(structure(list(covmat = covmat,
                          qvframe = frame,
                          dispersion = dispersion,
                          relerrs = relerrs,
                          factorname = factorname,
                          coef.indices = coef.indices,
                          modelcall = modelcall),
                     class="qv"))
}

indentPrint <- function(object, indent = 4, ...){
    zz <- ""
    tc <- textConnection("zz", "w", local = TRUE)
    sink(tc)
    try(print(object, ...))
    sink()
    close(tc)
    indent <- paste(rep(" ", indent), sep = "", collapse = "")
    cat(paste(indent, zz, sep = ""), sep = "\n")}

print.qv <- function(x, ...){
    print(x$qvframe)
}

summary.qv <- function(object, ...)
{
    if (!is.null(object$modelcall))
        cat("Model call: ",
            deparse(object$modelcall), "\n")
    if (!is.null(object$dispersion))
        cat("Dispersion: ", object$dispersion, "\n")
    if (!is.null(object$factorname))
        cat("Factor name: ",object$factorname,"\n")
    indentPrint(object$qvframe,...)
    if (!is.null(object$relerrs)){
        minErrSimple <- round(100*min(object$relerrs), 1)
        maxErrSimple <- round(100*max(object$relerrs), 1)
        errors<-worstErrors(object)
        minErrOverall<-round(100*errors[1], 1)
        maxErrOverall<-round(100*errors[2], 1)
        cat("Worst relative errors in SEs of simple contrasts (%): ",
            minErrSimple, maxErrSimple, "\n")
        cat("Worst relative errors over *all* contrasts (%): ",
            minErrOverall, maxErrOverall, "\n")
    }
}

plot.qv <- function(x,
                    intervalWidth = 2,
                    ylab = "estimate",
                    xlab = "",
                    ylim = NULL,
                    main = "Intervals based on quasi standard errors",
                    levelNames = NULL,
                    ...) {
    frame <- x$qvframe
    if (!is.null(levelNames)) {
        if (nrow(frame) != length(levelNames)) stop(
                                  "levelNames is not a vector of the right length"
                                               )
        row.names(frame) <- levelNames
    }
    if (is.null(frame$quasiSE))
        stop("Cannot plot, because there are no quasi standard errors")
    if (is.na(frame$estimate[1]))
        stop("No parameter estimates to plot")
    if (any(is.nan(frame$quasiSE)))
        stop(paste("No comparison intervals available,\n",
                   "since one of the quasi variances is negative.",
                   "  See ?qvcalc for more.",
                   sep = ""))
    faclevels <- factor(row.names(frame), levels = row.names(frame))
    xvalues <- seq(along = faclevels)
    est <- frame$estimate
    se <- frame$quasiSE
    tops <- est + (intervalWidth * se)
    tails <- est - (intervalWidth * se)
    range <- max(tops) - min(tails)
    if (is.null(ylim)) ylim <- c(min(tails) - range/10, max(tops) + range/10)
    plot(faclevels, frame$estimate, border = "transparent", ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, ...)
    points(frame$estimate, ...)
    segments(xvalues, tails, xvalues, tops)
    invisible(x)
}

