source("simulation.R")
# Simulation with culling -------------------------------------------------

np.culling.1    <- numeric(1000)
np.culling.2    <- numeric(1000)
np.culling.3    <- numeric(1000)
np.culling.5    <- numeric(1000)

for(j in 1:1000){
  
  id    <- sample(10000:100000, 1)
  f.bar <- f.matrix[, id]
  f     <- projection.matrix%*%f.bar
  f[d > 30] <- -20
  K     <- beta.vector.to.matrix(exp(f))
  
  
  new.outbreak              <- simulate.outbreak.with.culling(K, 3, 0.5, dist.mat, 1)
  np.culling.1[j]         <- new.outbreak$i.count
  
  new.outbreak              <- simulate.outbreak.with.culling(K, 3, 0.5, dist.mat, 2)
  np.culling.2[j]         <- new.outbreak$i.count
  
  new.outbreak              <- simulate.outbreak.with.culling(K, 3, 0.5, dist.mat, 3)
  np.culling.3[j]         <- new.outbreak$i.count
  
  new.outbreak              <- simulate.outbreak.with.culling(K, 3, 0.5, dist.mat, 5)
  np.culling.5[j]         <- new.outbreak$i.count
  
  
}


