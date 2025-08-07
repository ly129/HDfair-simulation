rm(list = ls())

seeds <- read.delim("seed_list.txt", sep = "\n", header = FALSE)

seeds <- as.vector(seeds[,1])

source("setting.R")
source("eval.R")

folder.name <- paste0("M",M,"A",A,"N",N/M,"p", p,"sdm", sd.th.m, "sda", sd.th.a, "sdy", sd.yi, "th", th.base[1], "x", p.nz) #, "wt", as.numeric(weighted))
folder.name
files <- paste0(folder.name, "/", seeds, "_res.rds")

mse.pred <- precision <- recall <- tp <- fp <- array(dim = c(A, 7, length(seeds)))
threshold <- 1e-3 * M * A
# threshold <- 0

for (i in 1:length(seeds)) {
  res <- readRDS(files[i])

  X.test <- res$X.test
  y.test <- res$y.test
  th.true <- res$th.true
  mse <- apply(res$estimates,
               MARGIN = 4,
               pred.error,
               ma = ma,
               xtest = X.test,
               ytest = y.test,
               simplify = TRUE)
  if (M == 3) {
    mse <- rbind(colSums(mse[1:3, ]), colSums(mse[4:6, ]), colSums(mse[7:9, ]))
  }
  mse.pred[,,i] <- mse/colSums(table(ma[,1], ma[,2]))

  # mse.esti <- apply(res$estimates,
  #                   MARGIN = 4,
  #                   esti.error,
  #                   th.true = res$th.true,
  #                   simplify = TRUE)

  for (j in 1:p) {
    res$estimates[j, , , "lam.min"] <- ifelse(rep(sum((res$estimates[j, , , "lam.min"])^2), M*A) < threshold,
                                              0,
                                              res$estimates[j, , , "lam.min"])
    res$estimates[j, , , "eta.min"] <- ifelse(rep(sum((res$estimates[j, , , "eta.min"])^2), M*A) < threshold,
                                              0,
                                              res$estimates[j, , , "eta.min"])
    for (m in 1:M) {
      for (a in 1:A) {
        res$estimates[j, a, m, "lasso s"] <- ifelse(res$estimates[j, a, m, "lasso s"]^2 < threshold/M/A,
                                                    0,
                                                    res$estimates[j, a, m, "lasso s"])
      }
    }

    for (m in 1:M) {
      res$estimates[j, , m, "lasso m"] <- ifelse(rep(sum((res$estimates[j, , m, "lasso m"])^2), A) < threshold/M,
                                                 0,
                                                 res$estimates[j, , m, "lasso m"])
    }

    for (a in 1:A) {
      res$estimates[j, a, , "lasso a"] <- ifelse(rep(sum((res$estimates[j, a, , "lasso a"])^2), M) < threshold/A,
                                                 0,
                                                 res$estimates[j, a, , "lasso a"])
    }

    res$estimates[j, , , "lasso p"] <- ifelse(rep(sum((res$estimates[j, , , "lasso p"])^2), M*A) < threshold,
                                              0,
                                              res$estimates[j, , , "lasso p"])
    res$estimates[j, , , "armul"] <- ifelse(rep(sum((res$estimates[j, , , "armul"])^2), M*A) < threshold,
                                            0,
                                            res$estimates[j, , , "armul"])
    # res$estimates[j, , , "RMTL"] <- ifelse(rep(sum((res$estimates[j, , , "RMTL"])^2), M*A) < threshold,
    #                                        0,
    #                                        res$estimates[j, , , "RMTL"])
  }

  precision.tmp <- apply(res$estimates,
                         MARGIN = 4,
                         prec,
                         simplify = TRUE,
                         th.true = th.true)
  if (M == 3) {
    precision.tmp <- rbind(colMeans(precision.tmp[1:3, ]), colMeans(precision.tmp[4:6, ]), colMeans(precision.tmp[7:9, ]))
  }
  precision[,,i] <- precision.tmp

  recall.tmp <- apply(res$estimates,
                      MARGIN = 4,
                      reca,
                      simplify = TRUE,
                      th.true = th.true)
  if (M == 3) {
    recall.tmp <- rbind(colMeans(recall.tmp[1:3, ]), colMeans(recall.tmp[4:6, ]), colMeans(recall.tmp[7:9, ]))
  }
  recall[,,i] <- recall.tmp

  tp.tmp <- apply(res$estimates,
                  MARGIN = 4,
                  truep,
                  simplify = TRUE,
                  th.true = th.true)

  if (M == 3) {
    tp.tmp <- rbind(colMeans(tp.tmp[1:3, ]), colMeans(tp.tmp[4:6, ]), colMeans(tp.tmp[7:9, ]))
  }
  tp[,,i] <- tp.tmp

  fp.tmp <-  apply(res$estimates,
                   MARGIN = 4,
                   falsep,
                   simplify = TRUE,
                   th.true = th.true)
  if (M == 3) {
    fp.tmp <- rbind(colMeans(fp.tmp[1:3, ]), colMeans(fp.tmp[4:6, ]), colMeans(fp.tmp[7:9, ]))
  }
  fp[,,i] <- fp.tmp
}


prediction.mse <- apply(mse.pred, 1:2, mean)
precision.mean <- apply(precision, 1:2, mean)
recall.mean <- apply(recall, 1:2, mean)
tp.mean <- apply(tp, 1:2, mean)
fp.mean <- apply(fp, 1:2, mean)

f1 <- 2*precision*recall/(precision+recall)
f1.mean <- apply(f1, 1:2, mean)

colnames(prediction.mse) <- colnames(precision.mean) <- colnames(recall.mean) <- colnames(f1.mean) <- colnames(tp.mean) <- colnames(fp.mean) <- dimnames(res$estimates)[[4]]


prediction.mse
precision.mean
recall.mean
f1.mean
tp.mean
fp.mean



out <- rbind(
  prediction.mse,
  precision.mean,
  recall.mean,
  f1.mean,
  tp.mean,
  fp.mean
)

out <- round(out, 2)

rownames(out) <- rep(c("prediction", "precision", "recall", "f1", "tp", "fp"), each = 3)

