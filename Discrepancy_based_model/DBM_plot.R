args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
beta.1 <- as.numeric(args[2])
beta.2 <- as.numeric(args[3])
INF <- as.numeric(args[4])
gamma.true <-  as.numeric(args[5])
ID <- as.numeric(args[6])
P <- 2
burn.in <- 2000

# Positions and Distances --------------------------------------------------
positions   <- read.csv(paste("positions", ID, ".dat", sep = ""), header = FALSE)
t           <- read.csv(paste("infecteds", ID, ".dat", sep = ""), header = FALSE)[, 1]
positions   <- matrix(positions[, 1], nrow = N, ncol = 2, byrow = TRUE)
dist.mat    <- dist(positions)                    #Construct Distance Matrix
d           <- sort(dist.mat)                     #Sort Distances
d0          <- read.csv(paste("distances0_", ID, ".txt", sep = ""), header = FALSE)[, 1]
d1          <- read.csv(paste("distances1_", ID, ".txt", sep = ""), header = FALSE)[, 1]
if(INF == 0){
  i <- read.csv(paste("i_times", ID, ".dat", sep = ""), header = FALSE)
}
beta.true <- as.numeric(read.table("beta_values.txt", header = FALSE))


# Import Matrices ---------------------------------------------------------

beta.mean0   <- read.csv(paste("mean0_", ID, ".txt", sep = ""), header = FALSE)[, 1]
beta.mean1   <- read.csv(paste("mean1_", ID, ".txt", sep = ""), header = FALSE)[, 1]
beta.var0    <- read.csv(paste("var0_", ID, ".txt", sep = ""), header = FALSE)[, 1]
beta.var1    <- read.csv(paste("var1_", ID, ".txt", sep = ""), header = FALSE)[, 1]
g            <- read.table(paste("gamma", ID, ".txt", sep = ""), header = FALSE)[, 1]
ell          <- read.table(paste("ell", ID, ".txt", sep = ""), header = FALSE)
beta.ci.0.upper <- beta.mean0 + 1.96*sqrt(beta.var0/4800)
beta.ci.0.lower <- beta.mean0 - 1.96*sqrt(beta.var0/4800)
if(INF == 0){
  i.sum    <- read.csv(paste("i_sum", ID, ".txt", sep = ""), header = FALSE)
  i.sum    <- i.sum[burn.in:(length(i.sum[, 1])), 1]
  i.sum.true <- sum(i[i < Inf])
}

beta.mean.1.approx <- approx(x = d1, y = beta.mean1, xout = d0)$y

# Plot --------------------------------------------------------------------
if(INF == 0){
  pdf(paste("Type_dist_plot_", ID, ".pdf", sep = ""))
  par(mfrow = c(2, 1))
  plot(d0, beta.mean0, type = 'l', main = "Infection Rate", xlab = "distance", ylab = expression(beta), ylim = c(0, max(beta.true[1], max(beta.mean0))))
  lines(d0, beta.true[1]*exp(-beta.1*d0), type = 'l', col = 'red')
  #polygon(c(d0, rev(d0)), c(beta.ci.0.lower, rev(beta.ci.0.upper)), col = rgb(1, 0, 0, 0.25), border = NA)
  plot(d1, beta.mean1, type = 'l', main = "Infection Rate", xlab = "distance", ylab = expression(beta), ylim = c(0, max(beta.true[2], max(beta.mean1))))
  lines(d1, beta.true[2]*exp(-beta.2*d1), type = 'l', col = 'red')
  legend("topright", bty = "n", legend = c("True Shape", "Posterior Mean"), lty = c(1, 1), col = c(2, 1))
  #rug(d.sparse)
  par(mfrow = c(2, 2))
  hist(g[-c(1:burn.in)], xlab = expression(gamma), freq = FALSE, main = "")
  abline(v = gamma.true, col = 2)
  plot(i.sum, type = 'l', xlab = "Iteration", ylab = "Sum of Infection Times")
  abline(h = i.sum.true, col = 'red')
  par(mfrow = c(2, 2))
  hist(ell$V1, xlab = "length scale 1", freq = FALSE, main = "")
  plot(ell$V1, type = 'l', xlab = "Iteration", ylab = "length scale 1")
  hist(ell$V2, xlab = "length scale 2", freq = FALSE, main = "")
  plot(ell$V2, type = 'l', xlab = "Iteration", ylab = "length scale 2")
  dev.off()
} else{
  pdf(paste("Type_dist_plot_", ID, ".pdf", sep = ""))
  par(mfrow = c(2, 1))
  plot(d0, beta.mean0, type = 'l', main = "Type 0 Infection Rate", xlab = "distance", ylab = expression(beta), ylim = c(0, max(beta.true[1], max(beta.mean0))))
  lines(d0, beta.true[1]*exp(-beta.1*d0), type = 'l', col = 'red')
  legend("topright", bty = "n", legend = c("True Shape", "Posterior Mean"), lty = c(1, 1), col = c(2, 1))
  plot(d1, beta.mean1, type = 'l', main = "Type 1 Infection Rate", xlab = "distance", ylab = expression(beta), ylim = c(0, max(beta.true[2], max(beta.mean1))))
  lines(d1, beta.true[2]*exp(-beta.2*d1), type = 'l', col = 'red')
  #abline(h = beta.true[2], col = 2)
  legend("topright", bty = "n", legend = c("True Shape", "Posterior Mean"), lty = c(1, 1), col = c(2, 1))
  hist(g[-c(1:burn.in)], xlab = expression(gamma), freq = FALSE, main = "")
  abline(v = gamma.true, col = 2)
  par(mfrow = c(2, 2))
  hist(ell$V1, xlab = "length scale 1", freq = FALSE, main = "")
  plot(ell$V1, type = 'l', xlab = "Iteration", ylab = "length scale 1")
  hist(ell$V2, xlab = "length scale 2", freq = FALSE, main = "")
  plot(ell$V2, type = 'l', xlab = "Iteration", ylab = "length scale 2")
  dev.off()
}
