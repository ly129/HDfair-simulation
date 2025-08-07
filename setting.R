M <- 3L
p <- 100L
# between site heterogeneity
sd.th.m <- 0.3
if (M == 1) sd.th.m <- 0
# between group heterogeneity
sd.th.a <- 0.2



A <- 3L
n.M <- rep(300, M)
N <- sum(n.M)
fixed.th <- TRUE
# weighted true: loss function is sum of per-task MSE
# weighted false: loss function is sum of per-task sum of squared errors
weighted <- FALSE





# used to generate true theta.
## true.theta = th.base + group heterogeneity + site heterogeneity
th.base <- rep(0.6, 10)
# th.base <- 1:10 * 0.1
p.nz <- length(th.base)

# proportions of groups
props <- matrix(c(0.7, 0.2, 0.1), nrow = M)
# props <- matrix(c(1/3, 1/3, 1/3), nrow = M)
# props <- matrix(c(0.6, 0.1, 0.1, 0.1, 0.1), nrow = M)
if (M > 1) {
  props <- matrix(props, nrow = M, ncol = A, byrow = TRUE)
  props[2, ] <- props[2, c(2, 1, 3)]
}
# group indicator
ma = matrix(0,N,2)
ma[,1] = rep(1:M, n.M)
for (m in 1:M) {
  ma[ma[, 1] == m, 2] <- rep(1:A, n.M[m] * props[m, ])
}

table(ma[, 1], ma[, 2])

# gaussian noise
sd.yi <- 1


