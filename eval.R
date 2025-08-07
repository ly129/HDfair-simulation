pred.error <- function(th, ma, xtest, ytest) {
  M <- dim(th)[3]
  A <- dim(th)[2]
  error <- matrix(nrow = M, ncol = A)
  for (m in 1:M) {
    for (a in 1:A) {
      thma <- th[, a, m]
      ma.id <- which(ma[, 1] == m & ma[, 2] == a)
      xma <- xtest[ma.id, ]
      yma <- ytest[ma.id]
      error[m, a] <- sum((yma - xma %*% thma)^2)
    }
  }
  return(error)
}

esti.error <- function(th, th.true) {
  apply((th - th.true)^2, 2:3, mean)
}

prec <- function(th, th.true) {
  p <- dim(th)[1]
  p.nz <- sum(th.true[,1,1] != 0)
  truth <- rep(FALSE, p)
  truth[1:p.nz] <- TRUE

  nz <- apply(th, 2:3, `!=`, 0)
  tp <- apply(apply(nz, 2:3, `&`, truth, simplify = TRUE), 2:3, sum)
  fp <- apply(apply(nz, 2:3, `&`, !truth, simplify = TRUE), 2:3, sum)

  precision <- tp/(tp+fp)
  return(precision)
}

reca <- function(th, th.true) {
  p <- dim(th)[1]
  p.nz <- sum(th.true[,1,1] != 0)
  truth <- rep(FALSE, p)
  truth[1:p.nz] <- TRUE

  nz <- apply(th, 2:3, `!=`, 0)
  tp <- apply(apply(nz, 2:3, `&`, truth, simplify = TRUE), 2:3, sum)
  fn <- apply(apply(!nz, 2:3, `&`, truth, simplify = TRUE), 2:3, sum)

  recall <- tp/(tp+fn)
  return(recall)
}


truep <- function(th, th.true) {
  p <- dim(th)[1]
  p.nz <- sum(th.true[,1,1] != 0)
  truth <- rep(FALSE, p)
  truth[1:p.nz] <- TRUE

  nz <- apply(th, 2:3, `!=`, 0)
  apply(apply(nz, 2:3, `&`, truth, simplify = TRUE), 2:3, sum)
}

falsep <- function(th, th.true) {
  p <- dim(th)[1]
  p.nz <- sum(th.true[,1,1] != 0)
  truth <- rep(FALSE, p)
  truth[1:p.nz] <- TRUE

  nz <- apply(th, 2:3, `!=`, 0)
  apply(apply(nz, 2:3, `&`, !truth, simplify = TRUE), 2:3, sum)
}