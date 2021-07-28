#GP Inference of the form beta_{ij} = f(d_{ij})
#Assumes fixed infection times

#Load functions
source("mcmc_functions.R")
source("simulation.R")


# Simulate outbreak -------------------------------------------------------
N <- 1000
coords <- cbind(runif(N), runif(N))
dist.mat <- as.matrix(dist(coords))

true.beta <- c(700/N, 0.7) #for Cpaper beta0 = 700/1000, beta1 = 0.7, alpha = 3, gamma = 0.5

set.seed(123456)
outbreak <- simulate.outbreak(true.beta, 3, 0.5, dist.mat)
outbreak$i.count


# Process data ------------------------------------------------------------
n <- outbreak$i.count
i <- outbreak$i.times
r <- outbreak$r.times

#Set susceptibles  to have i = r = Inf
i[is.na(i) == TRUE] <- Inf
r[is.na(r) == TRUE] <- Inf
infecteds <- which(i < Inf)



# Set up PPA --------------------------------------------------------------
d <- sort(dist(coords))
d.sparse <- c(seq(0, 10, 0.2), seq(15, 140, 5))
N.sparse <- length(d.sparse)
k <- sq.exp(d.sparse, d.sparse, 5, 4)
k.chol <- chol(k + 1e-10*diag(length(d.sparse)))
k.inv <- chol2inv(k.chol)

k.link <- sq.exp(d, d.sparse, 5, 4)
projection.matrix <- k.link%*%k.inv

#Matrix which goes from vector d to matrix of all pairwise distances
dist.mat.index <- matrix(0, N, N)
for(k in 1:N){
  for(j in 1:k){
    if(k != j){
      dist.mat.index[k, j] <- min(which(d == dist.mat[k, j]))
    }
  }
}
dist.mat.index <- dist.mat.index + t(dist.mat.index)


# Set up MCMC -------------------------------------------------------------

#Matrix of min(r_j, i_k) - min(i_k, i_j)
delta.min.matrix <- matrix(0, N, N)
for (k in infecteds)
  for (j in 1:N)
    delta.min.matrix[k, j] <- (min(r[k], i[j]) - min(i[k], i[j])) #difference of minimums

#For ech infected get the list of individuals infectious immediately before them
infected.lists <- vector(length = N, mode = "list")
for (j in infecteds) 
  infected.lists[[j]] <- intersect(which(r > i[j]), which(i < i[j]))

f.bar       <- rep(-5, N.sparse) #Set initial f
f           <- projection.matrix%*%f.bar
beta.matrix <- beta.vector.to.matrix(exp(f))
d.sum       <- double.sum.beta.update(beta.matrix, delta.min.matrix)
log.prod    <- log.product.beta.update(beta.matrix, infected.lists)
s.sum       <- sum(r[r < Inf] - i[i<Inf])
delta       <- 0.01 #tuning parameter
loglike      <- -d.sum + log.prod


# MCMC loop --------------------------------------------------------------------
n.iter   <- 100000
g        <- numeric(n.iter)
f.matrix <- matrix(NA, N.sparse, n.iter)
loglike.store <- numeric(n.iter)
pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
for(j in 1:n.iter){
  
  #Update f using underrelaxed
  f.bar.prop     <- sqrt(1-delta^2)*f.bar + delta*mvnorm.chol(N.sparse, k.chol)
  f.prop         <- projection.matrix%*%f.bar.prop
  beta.prop      <- beta.vector.to.matrix(exp(f.prop))
  d.sum.prop     <- double.sum.beta.update(beta.prop, delta.min.matrix)
  log.prod.prop  <- log.product.beta.update(beta.prop, infected.lists)
  loglike.prop   <- -d.sum.prop + log.prod.prop
  log.p.acc      <- loglike.prop - loglike
  
  if(log(runif(1)) < log.p.acc){ #Accept/reject step
    f.bar        <- f.bar.prop
    d.sum        <- d.sum.prop
    log.prod     <- log.prod.prop
    loglike      <- loglike.prop
  }
  
  #Update gamma 
  # g[j]           <- rgamma(1, 4*n + 1, 0.01 + s.sum)
  
  f.matrix[, j]  <- f.bar
  loglike.store[j] <- loglike
  setTxtProgressBar(pb, j) # update text progress bar after each iter
}



# Plot results ------------------------------------------------------------
f.bar.median <- apply(f.matrix[,-c(1:10000)], 1, median)
f.bar.upper <- apply(f.matrix[,-c(1:10000)], 1, quantile, 0.975)
f.bar.lower <- apply(f.matrix[,-c(1:10000)], 1, quantile, 0.025)
plot(d.sparse[d.sparse < 5], exp(f.bar.median[d.sparse < 5]), type = 'l', ylim = 1.1*c(0, max(exp(f.bar.upper[d.sparse < 5]), true.beta[1])), xlab = "distance (km)", ylab = "infection rate")
polygon(c(d.sparse, rev(d.sparse)), c(exp(f.bar.upper), rev(exp(f.bar.lower))), col = rgb(0, 0, 1, 0.25), border = NA)
lines(d.sparse, true.beta[1]*exp(-true.beta[2]*d.sparse), col = 2)


# 
# write.csv("~/OneDrive - The University of Nottingham/epidemics/posterior_samples/gp_model.csv", x = f.matrix, quote = FALSE, row.names = FALSE)
# 


# f.bar.mean <- rowMeans(f.matrix[,-c(1:1000)])
# f.mean    <- projection.matrix%*%f.bar.mean





# Simulate with culling ---------------------------------------------------




np.final.size <- numeric(1000)

for(i in 1:1000){
  id    <- sample(10000:100000, 1)
  f.bar <- f.matrix[, id]
  f     <- projection.matrix%*%f.bar
  f[d > 30] <- -20
  K     <- beta.vector.to.matrix(exp(f))
  
  new.outbreak              <- simulate.outbreak.beta.matrix(K, 3, 0.5)
  np.final.size[i]          <- new.outbreak$i.count

}


plot(density(exp.final.size, from = 0, to = 1000), xlab = "final size", main = "", ylim = c(0, 0.02), xlab = "final size")
lines(density(np.final.size, from = 0, to = 1000), lty = 2)
lines(density(logistic.final.size, from = 0, to = 1000), lty = 3)
lines(density(exp.final.size, from = 0, to = 1000), lty = 4)
legend("topright", legend = c("Exponential", "Nonparametric", "Logistic"), lty = c(1, 2, 3), bty = "n")



c("mean", mean(true.final.size), mean(logistic.final.size), mean(np.final.size))
c("median", median(true.final.size), median(logistic.final.size), median(np.final.size))

final.size.df <- data.frame("true"= true.final.size, "gp" = np.final.size, "logistic" = logistic.final.size, "exp" = exp.final.size)



png("final_size_distribution.png", width = 1280, height = 720)
par(mar=c(5,6,4,1)+.1)
plot(density(exp.final.size, from = 0, to = 1000), xlab = "Final size", main = "", ylim = c(0, 0.02), cex.lab =2.5, cex.axis = 2.5, lwd = 2)
lines(density(np.final.size, from = 0, to = 1000), lty = 2, lwd = 2)
lines(density(logistic.final.size, from = 0, to = 1000), lty = 3, lwd = 2)
legend("topleft", legend = c("Exponential", "GP", "Logistic"), lty = c(1, 2, 3), bty = "n", cex = 2.5, lwd = 2)
dev.off()
