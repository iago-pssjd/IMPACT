library(psych)

pedhazur <- structure(list(Group = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
2L), X = c(5L, 2L, 4L, 6L, 3L, 8L, 5L, 7L, 9L, 6L), Y = 1:10), .Names = c("Group", 
"X", "Y"), class = "data.frame", row.names = c(NA, -10L))
pedhazur
ped.stats <- statsBy(pedhazur,"Group")
print(ped.stats, short = FALSE)


(ssx_I <- sum((pedhazur$X[which(pedhazur$Group == 1)] - mean(pedhazur$X[which(pedhazur$Group == 1)]))^2))
(ssx_II <- sum((pedhazur$X[which(pedhazur$Group == 2)] - mean(pedhazur$X[which(pedhazur$Group == 2)]))^2))
(ssy_I <- sum((pedhazur$Y[which(pedhazur$Group == 1)] - mean(pedhazur$Y[which(pedhazur$Group == 1)]))^2))
(ssy_II <- sum((pedhazur$Y[which(pedhazur$Group == 2)] - mean(pedhazur$Y[which(pedhazur$Group == 2)]))^2))
(Sxy_I <- sum(((pedhazur$X - mean(pedhazur$X[which(pedhazur$Group == 1)]))*(pedhazur$Y - mean(pedhazur$Y[which(pedhazur$Group == 1)])))[which(pedhazur$Group == 1)]))
(Sxy_II <- sum(((pedhazur$X - mean(pedhazur$X[which(pedhazur$Group == 2)]))*(pedhazur$Y - mean(pedhazur$Y[which(pedhazur$Group == 2)])))[which(pedhazur$Group == 2)]))

(ssx_t <- sum((pedhazur$X - mean(pedhazur$X))^2))
var(pedhazur$X) == ssx_t/(nrow(pedhazur) - 1)

(ssy_t <- sum((pedhazur$Y - mean(pedhazur$Y))^2))
var(pedhazur$Y) == ssy_t/(nrow(pedhazur) - 1)





(Sxy_t <- sum((pedhazur$X - mean(pedhazur$X))*(pedhazur$Y - mean(pedhazur$Y))))
cov(pedhazur$X, pedhazur$Y) == Sxy_t/(nrow(pedhazur) - 1)

(b_I <- Sxy_I/ssx_I)
(b_II <- Sxy_II/ssx_II)
(b_t <- Sxy_t/ssx_t)


(rI <- Sxy_I/sqrt(ssx_I*ssy_I))
(rII <- Sxy_II/sqrt(ssx_II*ssy_II))
(rt <- Sxy_t/sqrt(ssx_t*ssy_t))
rt - cor(pedhazur$X, pedhazur$Y)

pedhazurG <- aggregate(cbind(X,Y) ~ Group, pedhazur, FUN = mean)

nj <- as.numeric(by(pedhazur, pedhazur$Group, FUN = nrow))

(ssx_b <- sum(nj*(pedhazurG$X - mean(pedhazur$X))^2))
(ssy_b <- sum(nj*(pedhazurG$Y - mean(pedhazur$Y))^2))
(Sxy_b <- sum(nj*(pedhazurG$X - mean(pedhazurG$X))*(pedhazurG$Y - mean(pedhazurG$Y))))
(rb <- Sxy_b/sqrt(ssx_b*ssy_b))

(ssx_w <- sum((pedhazur$X - ave(pedhazur$X, pedhazur$Group, FUN = mean))^2))
(ssy_w <- sum((pedhazur$Y - ave(pedhazur$Y, pedhazur$Group, FUN = mean))^2))
(Sxy_w <- sum((pedhazur$X - ave(pedhazur$X, pedhazur$Group, FUN = mean))*(pedhazur$Y - ave(pedhazur$Y, pedhazur$Group, FUN = mean))))
(rw <- Sxy_w/sqrt(ssx_w*ssy_w))

(rw2 <- sum(c(rI, rII)*nj/sum(nj)))



(etax <- ssx_b/ssx_t)
(etay <- ssy_b/ssy_t)
