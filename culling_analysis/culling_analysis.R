source("simulation.R")
# Simulation with culling -------------------------------------------------

np.culling   <- numeric(1000)

for(j in 1:1000){
  
  id    <- sample(10000:100000, 1) #Choose sample from posterior
  f.bar <- f.matrix[, id] 
  f     <- projection.matrix%*%f.bar #Get full function from MPA
  K     <- beta.vector.to.matrix(exp(f)) #turn into infection matrix
  
  #Simulate outbreak with culling
  new.outbreak              <- simulate.outbreak.with.culling(K, 3, 0.5, dist.mat, 3)
  np.culling[j]             <- new.outbreak$i.count
  
  
}


