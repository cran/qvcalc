#' @export
#'
#' @rdname qvcalc
#' 
qvcalc.lm <- function(object, factorname = NULL, coef.indices = NULL,
                      dispersion = NULL, ...)
{
    coef.indices.saved <- coef.indices
    if (is.null(factorname) && is.null(coef.indices)) {
          stop("arguments 'factorname' and 'coef.indices' are both NULL")
      }
    if (is.null(coef.indices)) {   ## try to use factorname
        term.index <- which(attr(terms(object),"term.labels") == factorname)
        modelmat <- model.matrix(object)
        has.coef <- colnames(modelmat) %in% names(coef(object))
        coef.indices <- which(attr(modelmat, "assign")[has.coef] == term.index)
        if (length(object$xlevels[[factorname]]) == length(coef.indices)){
            ## factor has no constraint applied, eg if no intercept in model
            contmat <- diag(length(coef.indices))}
        else
            contmat <- eval(call(object$contrasts[[factorname]],
                                 object$xlevels[[factorname]]))
        estimates <- contmat %*% coef(object)[coef.indices]
        covmat <- vcov(object, dispersion = dispersion)
        covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
        covmat <- contmat %*% covmat %*% t(contmat)
    }
    else {
        k <- length(coef.indices)
        refPos <- numeric(0)
        if (0 %in% coef.indices) { ## there's a reference level to include
            refPos <- which(coef.indices == 0)
            coef.indices <- coef.indices[-refPos]
        }
        covmat <- vcov(object, dispersion = dispersion)
        covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
        estimates <- coef(object)[coef.indices]
        if (length(refPos) == 1) {
            if (length(estimates) != k) estimates <- c(0, estimates)
            covmat <- cbind(0, rbind(0, covmat))
            names(estimates)[1] <- rownames(covmat)[1] <-
                colnames(covmat)[1] <- "(reference)"
            if (refPos != 1) {
                if (refPos == k){
                    perm <- c(2:k, 1)
                } else {
                    perm <- c(2:refPos, 1, (refPos + 1):k)
                }
                estimates <- estimates[perm]
                covmat <- covmat[perm, perm, drop = FALSE]
            }
        }
    }
    return(qvcalc.default(covmat,
                  factorname = factorname,
                  coef.indices = coef.indices.saved,
                  labels = NULL,
                  dispersion = dispersion,
                  estimates = estimates,
                  modelcall = object$call, ...)
           )
}
