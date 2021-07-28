#Functions to perform GP inference

double.sum <- function(i.times, r.times, beta.matrix, n, N) {
  #Calculate the integral of S_tI_t dt
  #Inputs: i.times - infection times, r.times - recovery_times, N - total population size, g - values of GP, d - sorted distances
  #dist.mat - matrix of distance, n <- length(i.times) #calculate number of infecteds
  
  # if(length(i.times) < N){ #if there are non-infecteds, their infection times are infitiny
  #   i.times[n+1:N] <- Inf
  # }
  result <- 0
  infecteds <- which(i < Inf)
  for (j in infecteds)
    for (k in 1:N)
      result <- result + beta.matrix[j, k]*(min(r.times[j], i.times[k]) - min(i.times[j], i.times[k])) #difference of minimums
  
  return(sum(result)) #sum values and return
}

double.sum.beta.update <- function(beta.matrix, delta.min.matrix) sum(beta.matrix*delta.min.matrix) #Update double sum if infection times are fixed


log.product <- function(n, i, r, beta.matrix){
  #Compute the log product term in the loglikelihood
  log.prod <- numeric(n)
  for (j in infecteds) {
    x <- numeric(1)
    x <- intersect(which(r > i[j]), which(i < i[j]))
    log.prod[j] <- sum(beta.matrix[x, j])
  }
  log.prod[which(log.prod == 0)] <- 1 #to preserve prods and logs
  log.prod[is.na(log.prod)] <- 1      #to preserve prods and logs
  
  return(sum(log(log.prod)))
}

log.product.beta.update <- function(beta.matrix, infection.list){
  #Update loglikelihood if infection times are fixed
  log.prod <- numeric(n)
  for (j in infecteds) {
    log.prod[j] <- sum(beta.matrix[infection.list[[j]], j])
  }
  log.prod[which(log.prod == 0)] <- 1 #to preserve prods and logs
  log.prod[is.na(log.prod)] <- 1      #to preserve prods and logs
  
  return(sum(log(log.prod)))
  
}


log.infectious.period.contributions <- function(alpha, g, n, s.sum) alpha*n*log(g) - g*s.sum

loglike.function <- function(d.sum, log.prod, log.pdf) d.sum + log.prod + log.pdf

sq.exp <- function(x, 
                   y, 
                   alpha, 
                   ell) {
  #Calculate the Covariance Matrix with the covariance function of squared exponential
  #INPUTS: x, y -- data points, alpha -- vertical scale parameter, ell -- horizontal scale
  ##OUTPUTS: covar -- length(x)*length(y) covariance matrix
  
  covar <- matrix(rep(0, length(x) * length(y)), nrow = length(x)) # Initialise Empty Covariance Matrix
  for (i in 1:length(x)) {
    for (j in 1:length(y))
      covar[i, j] <- alpha ^ 2 * exp(-(x[i] - y[j]) ^ 2 / ell ^ 2)
  }
  return(covar)
}


mvnorm.chol <- function(N, k.chol){
  #Multivariate Normal Sampler with Cholesky Input
  #Inputs: mu -- mean, chol -- cholesky decomposition of variance matrix (may need transposing)
  return(t(k.chol)%*%rnorm(N))  
}


beta.vector.to.matrix <- function(vec){
  #Turn vector into matrix
  beta.matrix <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:i){
      if(i != j)
        beta.matrix[i, j] <- vec[dist.mat.index[i, j]]
    }
  }
  beta.matrix <- beta.matrix + t(beta.matrix)
  return(beta.matrix)
}

two.cov.beta.vector.to.matrix <- function(vec, mat){
  #Turn vector into matrix
  beta.matrix <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:i){
      if(i != j)
        beta.matrix[i, j] <- vec[mat[i, j]]
    }
  }
  beta.matrix <- beta.matrix + t(beta.matrix)
  return(beta.matrix)
}



compute.dic <- function(loglike, mcmc.samples, burn.in, n, dist.mat, type){

  
  #loglike at posterior mean
  theta.hat         <- colMeans(mcmc.samples[-c(1:burn.in), ])
  if(type == "exp"){
    plain.K           <- exp(-theta.hat[2]*dist.mat)
  } else if(type == "logistic"){
      plain.K           <- 1/(1 + dist.mat^theta.hat[2])
  } else if(type == "homo"){
     plain.K            <- matrix(1, dim(dist.mat)[1], dim(dist.mat)[1])
  } else{
    stop("unknown type.")
  }

  plain.d.sum       <- double.sum.beta.update(plain.K, delta.min.matrix)
  d.sum             <- theta.hat[1]*plain.d.sum
  log.prod          <- (n-1)*log(theta.hat[1]) + log.product.beta.update(plain.K, infected.lists)
  loglike.at.mean   <- as.numeric(-d.sum + log.prod)
  
  pdic <-  2*(loglike.at.mean - mean(loglike[-c(1:burn.in)]))
  dic  <- -2*loglike.at.mean + 2*pdic
  
  return(dic)
  
}



