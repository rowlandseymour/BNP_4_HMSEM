args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
ID <- as.numeric(args[2])
S <- as.numeric(args[3])
x <- read.csv(paste("positions", ID, ".dat", sep = ""), header = FALSE)[, 1]
x <- matrix(x, nrow = N, byrow = TRUE)
t <- read.csv(paste("types", ID, ".dat", sep = ""), header = FALSE)[, 1]
dist.mat <- as.matrix(dist(x))
d0 <- dist.mat[which(t==0), ]
d0.sort <- unique(sort(d0))
d1 <- dist.mat[which(t==1), ]
d1.sort <- unique(sort(d1))
sparse.indexes.0 <- sample(2:(length(d0.sort)-1), S-2)
sparse.indexes.1 <- sample(2:(length(d1.sort)-1), S-2)
d0.sparse <- sort(c(min(d0), d0.sort[sparse.indexes.0], max(d0)))
d1.sparse <- sort(c(min(d1), d1.sort[sparse.indexes.1], max(d1)))
unique0 <-  length(d0.sort)
unique1 <- length(d1.sort)
lengths <- c(unique0, unique1)
index.mat.0 <- matrix(NA, nrow = length(which(t==0)), ncol = N)
index.mat.1 <- matrix(NA, nrow = length(which(t==1)), ncol = N)
for(i in 1:length(which(t==0))){
  for(k in 1:N){
      index.mat.0[i, k] <- which(d0.sort == d0[i, k])
  }
}
for(i in 1:length(which(t==1))){
  for(k in 1:N){
      index.mat.1[i, k] <- which(d1.sort == d1[i, k])
  }
}
index.mat.0 <- index.mat.0 - 1
index.mat.1 <- index.mat.1 - 1
index.mat.0 <- t(index.mat.0)
index.mat.1 <- t(index.mat.1)
write.table(lengths, paste("lengths_", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(d0.sort, paste("distances0_", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(d1.sort, paste("distances1_", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(d0.sparse, paste("distances0_sparse", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(d1.sparse, paste("distances1_sparse", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.numeric(t(index.mat.0)), paste("index_mat0_", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.numeric(t(index.mat.1)), paste("index_mat1_", ID, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
