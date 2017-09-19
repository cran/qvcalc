worstErrors <- function(qv.object)
{
    reducedForm <- function(covmat, qvmat){
        nlevels <- dim(covmat)[1]
 	firstRow <- covmat[1, ]
 	ones <- rep(1, nlevels)
 	J <- outer(ones, ones)
 	notzero <- 2:nlevels
 	r.covmat <- covmat + (firstRow[1]*J) -
            outer(firstRow, ones) -
            outer(ones, firstRow)
 	r.covmat <- r.covmat[notzero, notzero]
 	qv1 <- qvmat[1, 1]
 	r.qvmat <- (qvmat + qv1*J)[notzero, notzero]
 	list(r.covmat, r.qvmat)}
    covmat <- qv.object$covmat
    qvmat <- diag(qv.object$qvframe$quasiVar)
    r.form <- reducedForm(covmat, qvmat)
    r.covmat <- r.form[[1]]
    r.qvmat <- r.form[[2]]
    inverse.sqrt <- solve(chol(r.covmat))
    evalues <- eigen(t(inverse.sqrt) %*% r.qvmat %*% inverse.sqrt,
                     symmetric=TRUE)$values
    sqrt(c(min(evalues), max(evalues))) - 1
}
