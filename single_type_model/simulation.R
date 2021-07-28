################################################################
#       Heterogenously Mixing Stochastic Epidemic              #
#                                                              #
#                     16/11/16                                 #
#                                                              #
#                                                              #
################################################################



simulate.outbreak <- function(b, alpha, gamm, dist.mat){
  # data_simulation generates N points uniformaly on a square, and selects 1 to be infected.
  # It then chooses the next individual to become infected with transmission coefficient
  # b_ij = b0exp(-p(i, j)), where p(i, j) is the euclidean distance between i and j. 
  # The infectious period is from Gamma(alpha, gamma). So this is non-markovian.
  # Recovery times are generated at infections.
  # 
  # INPUTS: bet0 - transmission coefficient, 
  # alpha, gamma - infectious period coeffs
  # N - Population Size, plots - plotted output (default FALES)
  #
  # OUTPUTS: $i_times - infection times, $r_times - recovery times, $i_count - number of infecteds
  # Note: output times centred around min(r_times). 
  

  
  
  
  # Initial Vectors ------------------------------------------------------
  N <- dim(dist.mat)[1]
  s <- numeric(2*N) #number of susceptibles
  i <- numeric(2*N) # number of infecteds
  times <- data.frame(matrix(NA, nrow = N, ncol = 2)) #recovery and infection times
  time <- numeric(2*N) #all times
  R <- numeric(N) # recovery times to write over
  i_count <- 1 #count number of infecteds. We start with 1 inital infected
  r_count <- 0 #recovered count
  bet_vec <- numeric(N) #vector of beta values for each indiviudal
  index <- 1:N
  recovereds <- numeric(N) #vector of recovereds
  infecteds <- numeric(1) #vector of infecteds
  
  
  
  
  
  
  
  
  # Initial Conditions ------------------------------------------------------
  i0 <- sample(1:N, 1)
  i[1] <- 1 #1 infected Initially
  s[1] <- N - 1
  infecteds[1] <- i0 #1st infected is i0
  susceptibles <- c(1:N)[-i0]
  
  # Infectious Period Parameters --------------------------------------------
  
  times[i0, 1] <- 0 #start at time t = 0
  R[i0] <- rgamma(1, alpha, gamm) #generate first recovery time
  times[i0, 2] <- R[i0]
  
  
  
  beta.matrix <- b[1]*exp(-b[2]*dist.mat)
  
  
  
  
  # Simulation Loop ---------------------------------------------------------
  
  ninfects <- i[1] #while condition check
  j <- 1
  while(ninfects > 0){
    if(s[j]  > 0){ #if there are susceptibles check if infection
      
     bet_mat <- beta.matrix[infecteds, susceptibles]
      
      
      
      
      #Generate Infection Times
      if(length(infecteds) == 1 & length(susceptibles) == 1){
        times_prop <- rexp(1, sum(bet_mat))
        next_infect <- susceptibles
      } else if(length(infecteds) == 1){
        times_prop <- rexp(1, sum(bet_mat))
        next_infect <- sample(susceptibles, 1, prob = bet_mat/sum(bet_mat))
      } else if(length(susceptibles) == 1){
        times_prop <- rexp(length(infecteds), bet_mat)
        infector <- which(times_prop == min(times_prop))
        next_infect <- susceptibles
      } else{
        times_prop <- rexp(length(infecteds), rowSums(bet_mat))
        infector <- which(times_prop == min(times_prop))
        next_infect <- sample(susceptibles, 1, prob = as.numeric(bet_mat[infector, ]/sum(bet_mat[infector, ])))
      }
      time_new <- min(times_prop) + time[j]
      
      if(time_new < min(R[R>0])){ #if this time is before a recovery, and infection occurs
        times[next_infect, 1] <- time_new #set infection time
        times[next_infect, 2] <- time_new + rgamma(1, alpha, gamm) #generate new recovery time
        R[next_infect] <- times[next_infect, 2]
        s[j+1] <- s[j] - 1 #decrease suscpetibles by 1
        i[j+1] <- i[j] + 1 #increase infecteds by 1
        time[j+1] <- time_new #set new time
        i_count <- i_count + 1 #increase number of infecteds
        infecteds[length(infecteds)+1] <- next_infect
        susceptibles <- susceptibles[-which(susceptibles == next_infect)]
      } else { #If next infection is after recovery time, a recovery occurs
        s[j+1] <- s[j] #susceptibles remain the same
        i[j+1] <- i[j] - 1 #an infected enters the recovered class
        next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
        time[j+1] <- R[next_recover]
        R[next_recover] <- 0 #remove recovery time
        r_count <- r_count + 1
        recovereds[r_count] <- next_recover
        infecteds <- infecteds[-which(infecteds == next_recover)]
        
      }
    } else { #if s=0, then recovery must occur next
      s[j+1] <- s[j] #susceptibles remain the same
      i[j+1] <- i[j] - 1 #an infected enters the recovered class
      next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
      time[j+1] <- R[next_recover] #time is set to be the recovery time
      R[next_recover] <- 0 #remove recovery time
      r_count <- r_count + 1
      recovereds[r_count] <- next_recover
      infecteds <- infecteds[-which(infecteds == next_recover)]
    }
    ninfects <- i[j+1] #check while condition
    j <- j +1
  }
  
  times <- times - min(times[, 2], na.rm = TRUE)
  
  
  
  
  
  
  
  
  # Plots -------------------------------------------------------------------
  
  
  
  return(list('i.times' = times[, 1], 'r.times' = times[, 2], 'i.count' = i_count, 'time' = time))
}



simulate.outbreak.beta.matrix <- function(b, alpha, gamm){
  # data_simulation generates N points uniformaly on a square, and selects 1 to be infected.
  # It then chooses the next individual to become infected with transmission coefficient
  # b_ij = b0exp(-p(i, j)), where p(i, j) is the euclidean distance between i and j. 
  # The infectious period is from Gamma(alpha, gamma). So this is non-markovian.
  # Recovery times are generated at infections.
  # 
  # INPUTS: b - beta.matrix
  # alpha, gamma - infectious period coeffs
  # N - Population Size, plots - plotted output (default FALES)
  #
  # OUTPUTS: $i_times - infection times, $r_times - recovery times, $i_count - number of infecteds
  # Note: output times centred around min(r_times). 
  
  
  
  
  
  # Initial Vectors ------------------------------------------------------
  
  s <- numeric(2*N) #number of susceptibles
  i <- numeric(2*N) # number of infecteds
  times <- data.frame(matrix(NA, nrow = N, ncol = 2)) #recovery and infection times
  time <- numeric(2*N) #all times
  R <- numeric(N) # recovery times to write over
  i_count <- 1 #count number of infecteds. We start with 1 inital infected
  r_count <- 0 #recovered count
  bet_vec <- numeric(N) #vector of beta values for each indiviudal
  index <- 1:N
  recovereds <- numeric(N) #vector of recovereds
  infecteds <- numeric(1) #vector of infecteds
  
  
  
  
  
  
  
  
  # Initial Conditions ------------------------------------------------------
  i0 <- sample(1:N, 1)
  i[1] <- 1 #1 infected Initially
  s[1] <- N - 1
  infecteds[1] <- i0 #1st infected is i0
  susceptibles <- c(1:N)[-i0]
  
  # Infectious Period Parameters --------------------------------------------
  
  times[i0, 1] <- 0 #start at time t = 0
  R[i0] <- rgamma(1, alpha, gamm) #generate first recovery time
  times[i0, 2] <- R[i0]
  
  
  
  beta.matrix <- b
  
  
  
  
  # Simulation Loop ---------------------------------------------------------
  
  ninfects <- i[1] #while condition check
  j <- 1
  while(ninfects > 0){
    if(s[j]  > 0){ #if there are susceptibles check if infection
      
      bet_mat <- beta.matrix[infecteds, susceptibles]
      
      
      
      
      #Generate Infection Times
      if(length(infecteds) == 1 & length(susceptibles) == 1){
        times_prop <- rexp(1, sum(bet_mat))
        next_infect <- susceptibles
      } else if(length(infecteds) == 1){
        times_prop <- rexp(1, sum(bet_mat))
        next_infect <- sample(susceptibles, 1, prob = bet_mat/sum(bet_mat))
      } else if(length(susceptibles) == 1){
        times_prop <- rexp(length(infecteds), bet_mat)
        infector <- which(times_prop == min(times_prop))
        next_infect <- susceptibles
      } else{
        times_prop <- rexp(length(infecteds), rowSums(bet_mat))
        infector <- which(times_prop == min(times_prop))
        next_infect <- sample(susceptibles, 1, prob = as.numeric(bet_mat[infector, ]/sum(bet_mat[infector, ])))
      }
      time_new <- min(times_prop) + time[j]
      
      if(time_new < min(R[R>0])){ #if this time is before a recovery, and infection occurs
        times[next_infect, 1] <- time_new #set infection time
        times[next_infect, 2] <- time_new + rgamma(1, alpha, gamm) #generate new recovery time
        R[next_infect] <- times[next_infect, 2]
        s[j+1] <- s[j] - 1 #decrease suscpetibles by 1
        i[j+1] <- i[j] + 1 #increase infecteds by 1
        time[j+1] <- time_new #set new time
        i_count <- i_count + 1 #increase number of infecteds
        infecteds[length(infecteds)+1] <- next_infect
        susceptibles <- susceptibles[-which(susceptibles == next_infect)]
      } else { #If next infection is after recovery time, a recovery occurs
        s[j+1] <- s[j] #susceptibles remain the same
        i[j+1] <- i[j] - 1 #an infected enters the recovered class
        next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
        time[j+1] <- R[next_recover]
        R[next_recover] <- 0 #remove recovery time
        r_count <- r_count + 1
        recovereds[r_count] <- next_recover
        infecteds <- infecteds[-which(infecteds == next_recover)]
        
      }
    } else { #if s=0, then recovery must occur next
      s[j+1] <- s[j] #susceptibles remain the same
      i[j+1] <- i[j] - 1 #an infected enters the recovered class
      next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
      time[j+1] <- R[next_recover] #time is set to be the recovery time
      R[next_recover] <- 0 #remove recovery time
      r_count <- r_count + 1
      recovereds[r_count] <- next_recover
      infecteds <- infecteds[-which(infecteds == next_recover)]
    }
    ninfects <- i[j+1] #check while condition
    j <- j +1
  }
  
  times <- times - min(times[, 2], na.rm = TRUE)
  
  
  
  
  
  
  
  
  # Plots -------------------------------------------------------------------
  
  
  
  return(list('i.times' = times[, 1], 'r.times' = times[, 2], 'i.count' = i_count, 'time' = time))
}




simulate.outbreak.size.distance <- function(b, alpha, gamm, dist.mat, w){
  # data_simulation generates N points uniformaly on a square, and selects 1 to be infected.
  # It then chooses the next individual to become infected with transmission coefficient
  # b_ij = b0exp(-p(i, j)), where p(i, j) is the euclidean distance between i and j. 
  # The infectious period is from Gamma(alpha, gamma). So this is non-markovian.
  # Recovery times are generated at infections.
  # 
  # INPUTS: bet0 - transmission coefficient, 
  # alpha, gamma - infectious period coeffs
  # N - Population Size, plots - plotted output (default FALES)
  #
  # OUTPUTS: $i_times - infection times, $r_times - recovery times, $i_count - number of infecteds
  # Note: output times centred around min(r_times). 
  
  
  
  
  
  # Initial Vectors ------------------------------------------------------
  
  s <- numeric(2*N) #number of susceptibles
  i <- numeric(2*N) # number of infecteds
  times <- data.frame(matrix(NA, nrow = N, ncol = 2)) #recovery and infection times
  time <- numeric(2*N) #all times
  R <- numeric(N) # recovery times to write over
  i_count <- 1 #count number of infecteds. We start with 1 inital infected
  r_count <- 0 #recovered count
  bet_vec <- numeric(N) #vector of beta values for each indiviudal
  index <- 1:N
  recovereds <- numeric(N) #vector of recovereds
  infecteds <- numeric(1) #vector of infecteds
  
  
  
  
  
  
  
  
  # Initial Conditions ------------------------------------------------------
  i0 <- sample(1:N, 1)
  i[1] <- 1 #1 infected Initially
  s[1] <- N - 1
  infecteds[1] <- i0 #1st infected is i0
  susceptibles <- c(1:N)[-i0]
  
  # Infectious Period Parameters --------------------------------------------
  
  times[i0, 1] <- 0 #start at time t = 0
  R[i0] <- rgamma(1, alpha, gamm) #generate first recovery time
  times[i0, 2] <- R[i0]
  
  
  
  beta.matrix <- t(b[1]*exp(-b[2]*dist.mat)*sqrt(w/b[3]))
  
  
  
  
  # Simulation Loop ---------------------------------------------------------
  
  ninfects <- i[1] #while condition check
  j <- 1
  while(ninfects > 0){
    if(s[j]  > 0){ #if there are susceptibles check if infection
      
      bet_mat <- beta.matrix[infecteds, susceptibles]
      
      
      
      
      #Generate Infection Times
      if(length(infecteds) == 1 & length(susceptibles) == 1){
        times_prop <- rexp(1, sum(bet_mat))
        next_infect <- susceptibles
      } else if(length(infecteds) == 1){
        times_prop <- rexp(1, sum(bet_mat))
        next_infect <- sample(susceptibles, 1, prob = bet_mat/sum(bet_mat))
      } else if(length(susceptibles) == 1){
        times_prop <- rexp(length(infecteds), bet_mat)
        infector <- which(times_prop == min(times_prop))
        next_infect <- susceptibles
      } else{
        times_prop <- rexp(length(infecteds), rowSums(bet_mat))
        infector <- which(times_prop == min(times_prop))
        next_infect <- sample(susceptibles, 1, prob = as.numeric(bet_mat[infector, ]/sum(bet_mat[infector, ])))
      }
      time_new <- min(times_prop) + time[j]
      
      if(time_new < min(R[R>0])){ #if this time is before a recovery, and infection occurs
        times[next_infect, 1] <- time_new #set infection time
        times[next_infect, 2] <- time_new + rgamma(1, alpha, gamm) #generate new recovery time
        R[next_infect] <- times[next_infect, 2]
        s[j+1] <- s[j] - 1 #decrease suscpetibles by 1
        i[j+1] <- i[j] + 1 #increase infecteds by 1
        time[j+1] <- time_new #set new time
        i_count <- i_count + 1 #increase number of infecteds
        infecteds[length(infecteds)+1] <- next_infect
        susceptibles <- susceptibles[-which(susceptibles == next_infect)]
      } else { #If next infection is after recovery time, a recovery occurs
        s[j+1] <- s[j] #susceptibles remain the same
        i[j+1] <- i[j] - 1 #an infected enters the recovered class
        next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
        time[j+1] <- R[next_recover]
        R[next_recover] <- 0 #remove recovery time
        r_count <- r_count + 1
        recovereds[r_count] <- next_recover
        infecteds <- infecteds[-which(infecteds == next_recover)]
        
      }
    } else { #if s=0, then recovery must occur next
      s[j+1] <- s[j] #susceptibles remain the same
      i[j+1] <- i[j] - 1 #an infected enters the recovered class
      next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
      time[j+1] <- R[next_recover] #time is set to be the recovery time
      R[next_recover] <- 0 #remove recovery time
      r_count <- r_count + 1
      recovereds[r_count] <- next_recover
      infecteds <- infecteds[-which(infecteds == next_recover)]
    }
    ninfects <- i[j+1] #check while condition
    j <- j +1
  }
  
  times <- times - min(times[, 2], na.rm = TRUE)
  
  
  
  
   
  
  
  
  # Plots -------------------------------------------------------------------
  
  
  
  return(list('i.times' = times[, 1], 'r.times' = times[, 2], 'i.count' = i_count, 'time' = time))
}






simulate.outbreak.with.culling <- function(K, alpha, gamm, dist.mat, culling_radius){
  # data_simulation generates N points uniformaly on a square, and selects 1 to be infected.
  # It then chooses the next individual to become infected with transmission coefficient
  # b_ij = b0exp(-p(i, j)), where p(i, j) is the euclidean distance between i and j. 
  # The infectious period is from Gamma(alpha, gamma). So this is non-markovian.
  # Recovery times are generated at infections.
  # 
  # INPUTS: bet0 - transmission coefficient, 
  # alpha, gamma - infectious period coeffs
  # N - Population Size, plots - plotted output (default FALES)
  #
  # OUTPUTS: $i_times - infection times, $r_times - recovery times, $i_count - number of infecteds
  # Note: output times centred around min(r_times). 
  
  
  
  
  
  # Initial Vectors ------------------------------------------------------
  
  s <- numeric(2*N) #number of susceptibles
  i <- numeric(2*N) # number of infecteds
  times <- data.frame(matrix(NA, nrow = N, ncol = 2)) #recovery and infection times
  time <- numeric(2*N) #all times
  R <- numeric(N) # recovery times to write over
  i_count <- 1 #count number of infecteds. We start with 1 inital infected
  r_count <- 0 #recovered count
  bet_vec <- numeric(N) #vector of beta values for each indiviudal
  index <- 1:N
  recovereds <- numeric(N) #vector of recovereds
  infecteds <- numeric(1) #vector of infecteds
  
  
  
  
  
  
  
  
  # Initial Conditions ------------------------------------------------------
  i0 <- sample(1:N, 1)
  i[1] <- 1 #1 infected Initially
  s[1] <- N - 1
  infecteds[1] <- i0 #1st infected is i0
  susceptibles <- c(1:N)[-i0]
  
  # Infectious Period Parameters --------------------------------------------
  
  times[i0, 1] <- 0 #start at time t = 0
  R[i0] <- rgamma(1, alpha, gamm) #generate first recovery time
  times[i0, 2] <- R[i0]
  
  
  
  beta.matrix <- K
  
  
  
  
  # Simulation Loop ---------------------------------------------------------
  
  ninfects <- i[1] #while condition check
  j <- 1
  while(ninfects > 0){
    
    #Set abilities of the authorities
    if(i[j] < 30){
      multiplier <- 0
    } else if(i[j] < 60){
      multiplier <- 0.5
    } else{
      multiplier <- 1
    }
    
    
    if(s[j]  > 0){ #if there are susceptibles check if infection
      
      bet_mat <- beta.matrix[infecteds, susceptibles]
      

      
      
      #Generate Infection Times
      if(length(infecteds) == 1 & length(susceptibles) == 1){
        times_prop   <- rexp(1, sum(bet_mat))
        next_infect  <- susceptibles
      } else if(length(infecteds) == 1){
        times_prop   <- rexp(1, sum(bet_mat))
        next_infect  <- sample(susceptibles, 1, prob = bet_mat/sum(bet_mat))
      } else if(length(susceptibles) == 1){
        times_prop   <- rexp(length(infecteds), bet_mat)
        infector     <- which(times_prop == min(times_prop))
        next_infect  <- susceptibles
      } else{
        times_prop   <- rexp(length(infecteds), rowSums(bet_mat))
        infector     <- which(times_prop == min(times_prop))
        next_infect  <- sample(susceptibles, 1, prob = as.numeric(bet_mat[infector, ]/sum(bet_mat[infector, ])))
      }
      time_new <- min(times_prop) + time[j]
      
      if(time_new < min(R[R>0])){ #if this time is before a recovery, and infection occurs
        times[next_infect, 1] <- time_new #set infection time
        times[next_infect, 2] <- time_new + rgamma(1, alpha, gamm) #generate new recovery time
        R[next_infect] <- times[next_infect, 2]
        s[j+1] <- s[j] - 1 #decrease suscpetibles by 1
        i[j+1] <- i[j] + 1 #increase infecteds by 1
        time[j+1] <- time_new #set new time
        i_count <- i_count + 1 #increase number of infecteds
        infecteds[length(infecteds)+1] <- next_infect
        susceptibles <- susceptibles[-which(susceptibles == next_infect)]
      } else { #If next infection is after recovery time, a recovery occurs
        s[j+1] <- s[j] #susceptible remain the same
        i[j+1] <- i[j] - 1 #an infected enters the recovered class
        time_of_recovery <- min(R[R>0])
        next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
        time[j+1] <- R[next_recover]
        R[next_recover] <- 0 #remove recovery time
        r_count <- r_count + 1
        recovereds[r_count] <- next_recover
        infecteds <- infecteds[-which(infecteds == next_recover)]
        
        if(multiplier > 0 & culling_radius > 0){
          #Cull neighbours
          
          #Find neighbours to cull
          neighbours.to.cull <- which(dist.mat[next_recover, ] < multiplier*culling_radius)
          temp               <- length(neighbours.to.cull)
          neighbours.to.cull <- neighbours.to.cull[1:min(multiplier*100, temp)]
          
          number.culled         <- length(neighbours.to.cull) #count how many are culled
          i[j+1]                <- i[j+1] - number.culled #remove from infected numbers
          R[neighbours.to.cull] <- 0 #remove from list of possible removals
          r_count               <- r_count + number.culled
          if(number.culled > 0){
            infecteds             <-setdiff(infecteds, neighbours.to.cull) #remove from list of infecteds
           times[neighbours.to.cull, 2] <- time_of_recovery #overwrite removal times with new culling time
          }
        }
        
        
        
      }
    } else { #if s=0, then recovery must occur next
      s[j+1] <- s[j] #susceptibles remain the same
      i[j+1] <- i[j] - 1 #an infected enters the recovered class
      time_of_recovery <- min(R[R>0])
      next_recover <- index[R == min(R[R>0])] #individual with smallest recovery time
      time[j+1] <- R[next_recover] #time is set to be the recovery time
      R[next_recover] <- 0 #remove recovery time
      r_count <- r_count + 1
      recovereds[r_count] <- next_recover
      infecteds <- infecteds[-which(infecteds == next_recover)]
      
      
      if(multiplier > 0 & culling_radius > 0){
        #Cull neighbours
        
        #Find neighbours to cull
        neighbours.to.cull <- which(dist.mat[next_recover, ] < multiplier*culling_radius)
        temp               <- length(neighbours.to.cull)
        neighbours.to.cull <- neighbours.to.cull[1:min(multiplier*100, temp)]
        
        number.culled         <- length(neighbours.to.cull) #count how many are culled
        i[j+1]                <- i[j+1] - number.culled #remove from infected numbers
        R[neighbours.to.cull] <- 0 #remove from list of possible removals
        r_count               <- r_count + number.culled
        if(number.culled > 0){
          infecteds             <-setdiff(infecteds, neighbours.to.cull) #remove from list of infecteds
          times[neighbours.to.cull, 2] <- time_of_recovery #overwrite removal times with new culling time
        }
      }
      
      
    }
    ninfects <- i[j+1] #check while condition
    j <- j +1
  }
  
  times <- times - min(times[, 2], na.rm = TRUE)
  
  
  
  
  
  
  
  
  # Plots -------------------------------------------------------------------
  
  
  
  return(list('i.times' = times[, 1], 'r.times' = times[, 2], 'i.count' = i_count, 'time' = time))
}


