rm(list = ls())

library(HDfair)
library(glmnet)
library(MASS)
library(reticulate)

# download ARMUL from https://github.com/kw2934/ARMUL into the root folder

source_python("ARMUL-main/ARMUL.py")
np <- import("numpy")
test <- ARMUL("linear")
source("setting.R")

folder.name <- paste0("M",M,"A",A,"N",N/M,"p", p,"sdm", sd.th.m, "sda", sd.th.a,
                      "sdy", sd.yi, "th", th.base[1], "x", p.nz)
folder.name
if (!dir.exists(folder.name)) {
  dir.create(folder.name)
}

rho <- 1
adj <- 1
eps <- 1e-6
maxiter <- 1e4
verbose <- FALSE
lambda_length <- 20
eta_length <- 10
lambda_ratio <- ifelse(p == 100, 1e-2, 5e-2)
eta_ratio <- 1e-2
nfold <- 5


sim <- function(seed){
  message("Seed ", seed)
  set.seed(seed)

  if (fixed.th) {
    th.true <- array(0, dim = c(p,A,M))
    # th.true[1:p.nz, 1, 1] <- th.base
    # th.true[1:p.nz, 2, 1] <- th.base - sd.th.a
    # th.true[1:p.nz, 3, 1] <- th.base + sd.th.a
    for (m in seq(M)) {
      epsilon.m <- rep(0, p.nz)
      if (M > 1) {
        # epsilon.m <- runif(n = p.nz, min = 0, max = sd.th.m)
        epsilon.m <- rep(sd.th.m * (m - mean(1:M)), p.nz)
      }
      for (a in seq(A)) {
        # epsilon.a <- runif(n = p.nz, min = 0, max = sd.th.a)
        epsilon.a <- ((1:p.nz + a - 1) %% 3 - 1) * sd.th.a
        th.true[1:p.nz, a, m] <- th.base + epsilon.m + epsilon.a
      }
    }
  } else {
    ## th.true = th.base + group heterogeneity + site heterogeneity
    th.true <- array(0, dim = c(p,A,M))
    for (m in seq(M)) {
      epsilon.m <- rep(0, p.nz)
      if (M > 1) {
        # epsilon.m <- runif(n = p.nz, min = 0, max = sd.th.m)
        epsilon.m <- rnorm(n = p.nz, mean = 0, sd = sd.th.m)
      }
      for (a in seq(A)) {
        # epsilon.a <- runif(n = p.nz, min = 0, max = sd.th.a)
        epsilon.a <- rnorm(n = p.nz, mean = 0, sd = sd.th.a)     ########## shift the mean, fix sd
        th.true[1:p.nz, a, m] <- th.base + epsilon.m + epsilon.a
      }
    }
  }
  th.true[1:min(15, p.nz), ,]


  # generate X
  X <- X.test <- matrix(NA, N, p)
  y <- y.test <- numeric(N)

  # x, y and group (indicator) are all lists of length M
  for (m in 1:M) {
    for (a in 1:A) {
      nma <- n.M[m] * props[m, a]
      maids <- ma[, 1] == m & ma[, 2] == a

      # X[maids, ] <- mvrnorm(n = nma, mu = rep(0, p), Sigma = Sigma)
      X[maids, ] <- matrix(rnorm(nma * p, mean = 0, sd = 1), nrow = nma, ncol = p)
      y[maids] <- X[maids, ] %*% th.true[, a, m] + rnorm(nma, mean = 0, sd = sd.yi)

      # X.test[maids, ] <- mvrnorm(n = nma, mu = rep(0, p), Sigma = Sigma)
      X.test[maids, ] <- matrix(rnorm(nma * p, mean = 0, sd = 1), nrow = nma, ncol = p)
      y.test[maids] <- X.test[maids, ] %*% th.true[, a, m] + rnorm(nma, mean = 0, sd = sd.yi)
    }
  }

  cv_lambda <- HDfair_cv_lambda(
    X = X,
    y = y,
    ma = ma,
    lambda_length = lambda_length,
    lambda_ratio = lambda_ratio,
    nfold = nfold,
    eta = p * 1e6,
    rho = rho,
    adj=adj,
    eps=eps,
    maxiter=maxiter,
    verbose = verbose
  )

  pdf(file = paste0(folder.name, "/", seed, "sp_lambda.pdf"))
  par(mfrow = c(M, A))
  plot_sp_lambda(cv_lambda$sp, cv_lambda$lambda.min)
  dev.off()

  pdf(file = paste0(folder.name, "/", seed, "cv_lambda.pdf"))
  par(mfrow = c(1,1))
  plot_cv_lambda(cv_lambda)
  dev.off()

  # grid search in the nbhd of lambda.min
  bestlambdaid <- cv_lambda$index["lam.min"]

  lambda.range <- bestlambdaid + -2:2

  if (bestlambdaid > 18) lambda.range <- 16:20
  if (bestlambdaid < 4) lambda.range <- 2:6

  cv_etas <- vector(mode = "list", length = 5)

  for (i in 1:5) {
    lam.tmp <- cv_lambda$lambda[lambda.range[i]]

    cv_eta <- HDfair_cv_eta(
      X = X,
      y = y,
      ma = ma,
      lambda = lam.tmp,
      eta_length = eta_length,
      eta_ratio = eta_ratio,
      nfold = nfold,
      foldid = cv_lambda$foldid,
      rho = rho,
      adj=adj,
      eps=eps,
      maxiter=maxiter,
      verbose = verbose
    )

    pdf(file = paste0(folder.name, "/", seed, "cv_eta", i, ".pdf"))
    par(mfrow = c(1,1))
    plot_cv_eta(cv_eta)
    dev.off()

    # pdf(file = paste0(folder.name, "/", seed, "sp_eta", i, ".pdf"))
    # if (M == 1) par(mfrow = c(2*M, 5))
    # plot_sp_eta(cv_eta$sp, cv_eta$eta.min, 1:10)
    # dev.off()

    cv_etas[[i]] <- cv_eta
  }

  cvms <- lapply(cv_etas, function(cvcv) {min(cvcv$cvm)})

  cv_eta <- cv_etas[[which.min(cvms)]]

  pdf(file = paste0(folder.name, "/", seed, "cv_eta.pdf"))
  par(mfrow = c(1,1))
  plot_cv_eta(cv_eta)
  dev.off()

  pdf(file = paste0(folder.name, "/", seed, "sp_eta.pdf"))
  if (M == 1) par(mfrow = c(2*M, 5))
  plot_sp_eta(cv_eta$sp, cv_eta$eta.min, 1:10)
  dev.off()

  # ARMUL
  f.ma <- with(as.data.frame(ma), paste(V1, V2, sep = "_"))
  X.df <- as.data.frame(X)
  df_list <- split(X.df, f.ma)
  X_list<- lapply(df_list, function(dfi) as.matrix(dfi))
  y_list <- split(y, f.ma)

  X_np_list <- unname(lapply(X_list, function(mat) np$array(mat)))
  y_np_list <- unname(lapply(y_list, function(vec) np$array(vec)))

  data_train_py <- list(
    X_np_list,   # becomes data[0] in Python
    y_np_list    # becomes data[1] in Python
  )
  # alphas <- seq(0.1, 2, by = 0.1)   # 0.05, 0.10, ..., 0.50
  if (p == 100) {
    alphas <- exp(seq(log(0.1), log(100), length.out = 20))
  } else {
    alphas <- exp(seq(log(1e-3), log(1e0), length.out = 20))
  }

  # Build the Python‐side lambda list in one go
  lbd_list_py <- lapply(alphas, function(alpha) {
    np$array(alpha * sqrt(p / sapply(X_list, nrow)))
    # np$array(alpha * rep(1, 3))
  })
  test$cv(
    data_train_py,
    lbd_list  = lbd_list_py,  # your grid of λ’s
    model     = "vanilla",
    n_fold    = 5L,
    eta_global = 0.01,
    eta_local  = 0.01,
    T_global   = 500L,
    T_local    = 1L,
    intercept = FALSE
  )

  errors_cv <- py_to_r(test$errors_cv)

  pdf(file = paste0(folder.name, "/", seed, "cv_armul.pdf"))
  par(mfrow = c(1, 1))
  plot(x = alphas, y = errors_cv, log = "x")
  dev.off()
  # which.min(errors_cv)
  best_lbd_py <- test$lbd_list[[which.min(errors_cv)]]

  test$vanilla(
    data_train_py,
    best_lbd_py,
    eta_global = 0.01,
    eta_local  = 0.01,
    T_global   = 500L,
    T_local    = 1L,
    intercept  = FALSE
  )

  if (M == 1) {
    th.armul <- t(test$models$vanilla[,,1])
  } else {
    th.armul <- array(t(test$models$vanilla[,,1]), dim = c(p,A,M))
  }

  pdf(file = paste0(folder.name, "/", seed, "lasso.s.pdf"))
  th.lasso.s <- th.lasso.s.1se <- array(NA, dim = c(p,A,M))
  if (M == 1) par(mfrow = c(M*A, 2))
  for (m in 1:M) {
    for (a in 1:A) {
      maid <- ma[, 1] == m & ma[, 2] == a
      lasso.s <- cv.glmnet(x = X[maid, ],
                           y = y[maid],
                           nfolds = 5,
                           foldid = cv_lambda$foldid[maid],
                           standardize = FALSE,
                           intercept = FALSE)
      # cvs.lasso.s[[r]][[(m-1) * A + a]] <- lasso.s
      th.lasso.s[, a, m] <- as.numeric(coef(lasso.s, s = "lambda.min"))[-1]
      th.lasso.s.1se[, a, m] <- as.numeric(coef(lasso.s, s = "lambda.1se"))[-1]
      plot(lasso.s)
      plot(lasso.s$glmnet.fit, xvar = "lambda")
      abline(v = log(lasso.s$lambda.min))
    }
  }
  dev.off()
  # if (M == 1) th.lasso.s <- th.lasso.s[, , 1]

  th.lasso.p <- th.lasso.p.1se <- array(NA, dim = c(p,A,M))
  lasso.p <- cv.glmnet(x = X,
                       y = y,
                       nfolds = 5,
                       foldid = cv_lambda$foldid,
                       standardize = FALSE,
                       intercept = FALSE)
  for (m in 1:M) {
    for (a in 1:A) {
      th.lasso.p[, a, m] <- as.numeric(coef(lasso.p, s = "lambda.min"))[-1]
      th.lasso.p.1se[, a, m] <- as.numeric(coef(lasso.p, s = "lambda.1se"))[-1]
    }
  }
  pdf(file = paste0(folder.name, "/", seed, "lasso.p.pdf"))
  par(mfrow = c(1, 2))
  plot(lasso.p)
  plot(lasso.p$glmnet.fit, xvar = "lambda")
  abline(v = log(lasso.p$lambda.min))
  dev.off()

  pdf(file = paste0(folder.name, "/", seed, "lasso.m.pdf"))
  th.lasso.m <- th.lasso.m.1se <- array(NA, dim = c(p,A,M))
  par(mfrow = c(M, 2))
  for (m in 1:M) {
    mid <- ma[, 1] == m
    lasso.m <- cv.glmnet(x = X[mid, ],
                         y = y[mid],
                         nfolds = 5,
                         foldid = cv_lambda$foldid[mid],
                         standardize = FALSE,
                         intercept = FALSE)
    plot(lasso.m)
    plot(lasso.m$glmnet.fit, xvar = "lambda")
    abline(v = log(lasso.m$lambda.min))
    for (a in 1:A) {
      th.lasso.m[, a, m] <- as.numeric(coef(lasso.m, s = "lambda.min"))[-1]
      th.lasso.m.1se[, a, m] <- as.numeric(coef(lasso.m, s = "lambda.1se"))[-1]
    }
  }
  dev.off()

  pdf(file = paste0(folder.name, "/", seed, "lasso.a.pdf"))
  th.lasso.a <- th.lasso.a.1se <- array(NA, dim = c(p,A,M))
  par(mfrow = c(A, 2))
  for (a in 1:A) {
    aid <- ma[, 2] == a
    lasso.a <- cv.glmnet(x = X[aid, ],
                         y = y[aid],
                         nfolds = 5,
                         foldid = cv_lambda$foldid[aid],
                         standardize = FALSE,
                         intercept = FALSE)
    plot(lasso.a)
    plot(lasso.a$glmnet.fit, xvar = "lambda")
    for (m in 1:M) {
      th.lasso.a[, a, m] <- as.numeric(coef(lasso.a, s = "lambda.min"))[-1]
      th.lasso.a.1se[, a, m] <- as.numeric(coef(lasso.a, s = "lambda.1se"))[-1]
    }
  }
  dev.off()


  # # RMTL
  # ma_pairs <- paste0("m", ma[,1], "a", ma[,2])
  # Xlist <- lapply(split(as.data.frame(X), ma_pairs), as.matrix)
  # ylist <- lapply(split(y, ma_pairs), matrix, ncol = 1)
  # cvMTL_fit <- RMTL::cvMTL(
  #   X = Xlist,
  #   Y = ylist,
  #   type = "Regression",
  #   Regularization = "L21",
  #   Lam1_seq = exp(seq(log(1e1), log(1e-2), length.out = 20))
  # )
  #
  # fitMTL <- RMTL::MTL(
  #   X = Xlist,
  #   Y = ylist,
  #   type = "Regression",
  #   Regularization = "L21",
  #   Lam1 = cvMTL_fit$Lam1.min
  # )
  # th.rmtl <- fitMTL$W
  #
  # pdf(file = paste0(folder.name, "/", seed, "RMTL.pdf"))
  # par(mfrow = c(1, 1))
  # plot(cvMTL_fit)
  # dev.off()

  estimates <- array(dim = c(p,A,M,7))
  dimnames(estimates) <- list(paste0("V", 1:p),
                              paste0("A", 1:A),
                              paste0("M", 1:M),
                              c("lam.min",
                                "eta.min",
                                # "lam.1se", "eta.1d",
                                "lasso s", #"lasso s 1se" ,
                                "lasso m", #"lasso m 1se",
                                "lasso a", #"lasso a 1se",
                                "lasso p", #"lasso p 1se",
                                "armul"))
  estimates[,,,1] <- cv_lambda$sp$estimates[,,,cv_lambda$index["lam.min"]]
  estimates[,,,2] <- cv_eta$sp$estimates[,,,cv_eta$index["eta.min"]]
  # estimates[,,,3] <- cv_lambda$sp$estimates[,,,cv_lambda$index["lam.1se"]]
  # estimates[,,,4] <- cv_eta_1d$sp$estimates[,,,cv_eta_1d$index["eta.min"]]
  estimates[,,,3] <- th.lasso.s
  # estimates[,,,4] <- th.lasso.s.1se
  estimates[,,,4] <- th.lasso.m
  # estimates[,,,6] <- th.lasso.m.1se
  estimates[,,,5] <- th.lasso.a
  # estimates[,,,8] <- th.lasso.a.1se
  estimates[,,,6] <- th.lasso.p
  # estimates[,,,10] <- th.lasso.p.1se
  estimates[,,,7] <- th.armul
  # estimates[,,,8] <- th.rmtl
  iters.lam <- cv_lambda$sp$iterations
  iters.eta <- cv_eta$sp$iterations
  # iters.lam2 <- cv_lambda2$sp$iterations
  # iters.eta2 <- cv_eta2$sp$iterations

  res <- list(estimates = estimates,
              iterations.lam = iters.lam,
              iterations.eta = iters.eta,
              # iterations.lam2 = iters.lam2,
              # iterations.eta2 = iters.eta2,
              X = X,
              y = y,
              X.test = X.test,
              y.test = y.test,
              th.true = th.true)
}






args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a seed.")
seed <- as.integer(args[1])

cat("Running with seed:", seed, "\n")

res <- sim(seed)


saveRDS(res, file = paste0(folder.name, "/", seed, "_res.rds"))

# saveRDS(rbind(mse.pred, mse.esti, precision.pred, recall.pred, tp.pred, fp.pred),
#         file = paste0(folder.name, "/", seed, "_sum.rds"))
