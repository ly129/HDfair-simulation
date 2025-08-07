# Code for Single-source dense scenario in the manuscript
rm(list = ls())

library(HDfair)
library(reticulate)

source_python("ARMUL-main/ARMUL.py")
np <- import("numpy")
test <- ARMUL("linear")


######################################################################
# Modify number of replications and/or heterogeneity level as needed.

# Number of replications
n.rep <- 10

# Heterogeneity level # No, Low, Medium, High
heterogeneity <- "No"

######################################################################







set.seed(20250616)

A <- 3L
M <- 1L
p.nz <- 10L
p <- 10L
n <- 300L

# generate true theta
th <- matrix(0, nrow = p, ncol = A)
if (heterogeneity == "No") {
  # no heterogeneity
  th[, 1] <- 0.6
  th[, 2] <- 0.6
  th[, 3] <- 0.6
} else if (heterogeneity == "Low") {
  # low heterogeneity
  th[, 1] <- rep(c(0.5, 0.6, 0.7), length.out = p.nz)
  th[, 2] <- rep(c(0.6, 0.7, 0.5), length.out = p.nz)
  th[, 3] <- rep(c(0.7, 0.5, 0.6), length.out = p.nz)
} else if (heterogeneity == "Medium") {
  # medium heterogeneity
  th[, 1] <- rep(c(0.4, 0.6, 0.8), length.out = p.nz) # 0.5
  th[, 2] <- rep(c(0.6, 0.8, 0.4), length.out = p.nz) # 0.2
  th[, 3] <- rep(c(0.8, 0.4, 0.6), length.out = p.nz) # 0.8
} else if (heterogeneity == "High") {
  # high heterogeneity
  th[, 1] <- rep(c(0.3, 0.6, 0.9), length.out = p.nz) # 0.5
  th[, 2] <- rep(c(0.6, 0.9, 0.3), length.out = p.nz) # 0.2
  th[, 3] <- rep(c(0.9, 0.3, 0.6), length.out = p.nz) # 0.8
}
# if (heterogeneity == "No") {
#   # no heterogeneity
#   th[1:p.nz, 1] <- 0.6
#   th[1:p.nz, 2] <- 0.6
#   th[1:p.nz, 3] <- 0.6
# } else if (heterogeneity == "Low") {
#   # low heterogeneity
#   th[1:p.nz, 1] <- 0.6
#   th[1:p.nz, 2] <- 0.7
#   th[1:p.nz, 3] <- 0.5
# } else if (heterogeneity == "Medium") {
#   # medium heterogeneity
#   th[1:p.nz, 1] <- 0.6
#   th[1:p.nz, 2] <- 0.8
#   th[1:p.nz, 3] <- 0.4
# } else if (heterogeneity == "High") {
#   # high heterogeneity
#   th[1:p.nz, 1] <- 0.6
#   th[1:p.nz, 2] <- 0.9
#   th[1:p.nz, 3] <- 0.3
# }
th


# gaussian noise
sd.yi <- 1

# generate groups
pts <- c(0.7, 0.2, 0.1)

ma <- matrix(nrow = n, ncol = 2)
ma[, 1] <- 1
ma[, 2] <- rep(1:3, n * pts)


n.etas <- 30
etas <- exp(seq(log(1e0), log(1e-4), length.out = n.etas))
sp.etas <- lms <- armuls <- vector(mode = "list", length = n.rep)

alphas <- exp(seq(log(0.001), log(10), length.out = n.etas))

## store data
train.list <- test.list <- vector(mode = "list", length = n.rep)

# lam <- 10^(seq(0.5, -3.5, -0.1))
lam <- 0
rho <- 1
maxit <- 1e4
eps <- 1e-6

for (r in 1:n.rep) {
  message("Replication ", r)
  # generate X
  X <- matrix(rnorm(n = p * n, mean = 0, sd = 1), nrow = n)

  X.test <- matrix(rnorm(n = p * n, mean = 0, sd = 1), nrow = n)

  y <- y.test <- numeric(n)

  for (a in 1:A) {
    ida <- ma[, 2] == a
    Xa <- X[ida, ]
    Xa.test <- X.test[ida, ]
    tha <- th[, a]
    y[ida] <- c(Xa %*% tha) + rnorm(sum(ida), mean = 0, sd = sd.yi)
    y.test[ida] <- c(Xa.test %*% tha) + rnorm(sum(ida), mean = 0, sd = sd.yi)
  }

  train.list[[r]]$x <- X; train.list[[r]]$y <- y
  test.list[[r]]$x <- X.test; test.list[[r]]$y <- y.test

  # linear model
  fit.lm <- lm(y ~ X - 1)
  lms[[r]] <- fit.lm

  # th.init <- array(dim = c(p, M, A))
  # for (m in 1:M) {
  #   th.init[, m ,] <- fit.single.m$coefficients
  # }

  fit <- HDfair_sp_eta(
    X = X,
    y = y,
    ma = ma,
    lambda = lam,
    eta_length = n.etas,
    eta_seq = etas,
    rho = rho,
    weighted = FALSE,
    maxiter = maxit,
    verbose = FALSE,
    eps = eps
  )

  # par(mfrow = c(2, 5))
  # for (j in 1:p.nz) {
  #   matplot(x = etas, y = t(fit$estimates[j, , 1, ]), type = "l", log = "x")
  #   abline(h = fit.lm$coefficients[j])
  # }
  sp.etas[[r]] <- fit



  # ARMUL -------------------------------------------------------------------
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

  lbd_list_py <- lapply(alphas, function(alpha) {
    np$array(alpha * sqrt(p / sapply(X_list, nrow)))
    # np$array(alpha * rep(1, 3))
  })

  theta_path  <- vector("list", length(n.etas))

  # 5. loop to fit vanilla ARMUL at each λ
  for (i in 1:n.etas) {
    test$vanilla(
      data_train_py,
      lbd_list_py[[i]],
      eta_global = 0.01,
      eta_local  = 0.01,
      T_global   = as.integer(500),
      T_local    = as.integer(1),
      intercept  = FALSE
    )
    # extract and convert back to R
    theta_path[[i]] <- py_to_r(test$models$vanilla)[,,1]         # d × m matrix
  }

  th.armul <- simplify2array(theta_path)

  # for (j in 1:p.nz) {
  #   matplot(x = alphas, y = t(th.armul[, j, ]), type = "l", log = "x")
  #   abline(h = fit.lm$coefficients[j])
  # }

  armuls[[r]] <- th.armul
}

# results
## MSE
mse.hdfair <- mse.armul <- array(dim = c(A, n.rep, n.etas))
mse.lm <- matrix(nrow = n.rep, ncol = A)

for (r in 1:n.rep) {
  X.test <- test.list[[r]]$x
  y.test <- test.list[[r]]$y

  for (a in 1:A) {
    ida <- ma[, 2] == a
    Xa.test <- X.test[ida, ]
    ya.test <- y.test[ida]

    mse.lm[r, a] <- mean((Xa.test %*% lms[[r]]$coefficients - ya.test)^2)

    for (nn in 1:n.etas) {
      th.eta.nn <- sp.etas[[r]]$estimates[, a, 1, nn]
      mse.hdfair[a, r, nn] <- mean((Xa.test %*% th.eta.nn - ya.test)^2)

      mse.armul[a, r, nn] <- mean((Xa.test %*% armuls[[r]][a, , nn] - ya.test)^2)
    }
  }
}

mse.hdfair.mean <- apply(mse.hdfair, MARGIN = c(1, 3), FUN = mean, simplify = FALSE)
mse.armul.mean <- apply(mse.armul, MARGIN = c(1, 3), FUN = mean, simplify = FALSE)
mse.lm.mean <- colMeans(mse.lm)

par(mfrow = c(1, 1))
matplot(x = etas,
        y = t(mse.hdfair.mean),
        type = "l",
        log = "x",
        xlab = expression(eta),
        ylab = "Mean Squared Error",
        lwd = 2,
        col = 2:4,
        lty = rep(1, A),
        main = paste0(heterogeneity, " Heterogeneity"),
        # sub = paste0("HDfair, sd.y = ", sd.yi),
        # xlim = c(1e-5, 1e1),
        # ylim = c(0.9, 2.4),
)
# abline(h = mse.lm.mean, col = 2:4, lwd = 2)

legend("topright",
       # x = 1e-4,
       # y = 3,
       legend = c("Group 1 (Majority)", "Group 2", "Group 3"),
       lwd = 2,
       lty = 1,
       col = 2:4,
       bty = "n")

matplot(x = alphas,
        y = t(mse.armul.mean),
        type = "l",
        log = "x",
        xlab = expression(alpha),
        ylab = "Mean Squared Error",
        lwd = 2,
        col = 2:4,
        lty = rep(1, A),
        main = paste0(heterogeneity, " Heterogeneity"),
        sub = "ARMUL"
        # xlim = c(1e-5, 1e1),
        # ylim = c(1, 4),
)

abline(h = mse.lm.mean, col = 2:4, lwd = 2)

legend("topright",
       # x = 1e-4,
       # y = 3,
       legend = c("Group 1 (Majority)", "Group 2", "Group 3"),
       lwd = 2,
       lty = 1,
       col = 2:4,
       bty = "n")

file.name <- paste0("~/Library/CloudStorage/Box-Box/FairReg/sim_redo/", heterogeneity, "Heterogeneity_sdy", sd.yi, ".Rdata")
# save.image(file.name)
