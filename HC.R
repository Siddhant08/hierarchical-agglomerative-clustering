#Task-1
#read labels
y <- read.table('/home/cim/pgt/mhac016/CS5100_labs/label.txt',strip.white = TRUE)
#convert list to double
y <- as.numeric(unlist(y))

#read file
data=read.table('/home/cim/pgt/mhac016/CS5100_labs/nci.data.txt', header = FALSE, sep = "", dec = ".")
#transpose data as the file has 64 observations and 6830 features
data=t(data)
#normalize data
data <- scale(data)

#function for HAC
hac <- function(data,meth)
{ 
  #decide method for clustering
  switch(meth,'single'={meth='min'},'average'={meth='mean'},'complete'={meth='max'})
  
  #get distance values
  dist_matrix <- as.matrix(dist(data))
  #set 0 as inf in distance matrix
  diag(dist_matrix) <- Inf
  
  N <- dim(data)[1]
  
  #set rownames and colnames
  colnames(dist_matrix) <- -(1:N)
  rownames(dist_matrix) <- -(1:N)
  
  diag(dist_matrix) <- Inf
  merge <- matrix(0, N-1, 2)
  height <- vector(length = N-1)
  
  for (m in 1:(N-1)) { 
    cols <- colnames(dist_matrix)
    #get index of the minimum distance
    #min_dist_index <- which(new_mat==min(new_mat),arr.ind=TRUE)
    min_dist_index <- which(dist_matrix == min(dist_matrix), arr.ind = TRUE)[1,,drop=FALSE]
    height[m] <- min(dist_matrix) # The height is the value of the pair with the minimum distance
    
    # The row and column position of the minimum pair is stored as sequence m in the merge object
    merge[m,] <- as.numeric(cols[min_dist_index])
    #get clusters
    clusters <- c(min_dist_index, which(cols %in% cols[min_dist_index[1, cols[min_dist_index] > 0]]))
    colnames(dist_matrix)[clusters] <- m
    
    #merge closest pairs
    sl <- apply(dist_matrix[min_dist_index,], 2, meth)
    dist_matrix[min(min_dist_index),] <- sl
    dist_matrix[,min(min_dist_index)] <- sl
    
    #set minimum distance pair to inf so they are not used again
    diag(dist_matrix) <- Inf
    
    #in the final iteration, all elements of the distance matrix(dist_matrix) will be set to Inf
    dist_matrix[max(min_dist_index),] <- Inf
    dist_matrix[,max(min_dist_index)] <- Inf
  }
  
  #get clusters
  s <- list("merged_datapoints"=merge,"dendrogram_height"=height,"clusters"=clusters)
  return(cutree(s,14))
}

#Task-2
#call the above function to implement the clustering algorithm on the data
hac_single <- hac(data,'single')
hac_avg <- hac(data,'average')
hac_complete <- hac(data,'complete')

#create plots for all linkages
plot(hac_complete,main="Complete Linkage",xlab="",sub="",cex=1,pch=20,col=hac_complete+1)
plot(hac_avg,main="Average Linkage",xlab="",sub="",cex=1,pch=20,col=hac_avg+1)
plot(hac_single,main="Single Linkage",xlab="",sub="",cex=1,pch=20,col=hac_single+1)

#Task-3
#check performance of the clustering algorithms
#need to install mclust package first. run the below statements in the console first.
#install.packages("mclust");
library(mclust)
cat("The similarity of clusters with class labels for single linkage: ",adjustedRandIndex(y,hac_single))
cat("The similarity of clusters with class labels for average linkage: ",adjustedRandIndex(y,hac_avg))
cat("The similarity of clusters with class labels for complete linkage: ",adjustedRandIndex(y,hac_complete))

#Task-4
#apply k-means
#perform k means clustering with k=14
set.seed(2)
km.out <- kmeans(data, 14, nstart=20)
#cluster assignments of the samples
km.out$cluster
#model performance
km.out$tot.withinss
cat("The similarity of clusters with class labels for k-means: ",adjustedRandIndex(y,km.out$cluster))

#plot the data points with different colours for the clusters
plot(data, col=(km.out$cluster+1), main="K-Means Clustering Results with K=14", xlab="", ylab="", pch=20, cex=2)

#trying k-means with different values of k
t <- Sys.time()
best_score <- 0
best_k_val <- 0
for (i in 1:14)
{
  km.out <- kmeans(data, i, nstart=20)
  #cluster assignments of the samples
  km.out$cluster
  #model performance
  km.out$tot.withinss
  #cat("The similarity of clusters with class labels for complete linkage: ",adjustedRandIndex(y,km.out$cluster))
  if (best_score<adjustedRandIndex(y,km.out$cluster))
  {
    best_score <- adjustedRandIndex(y,km.out$cluster)
    best_k_val <- i
  }
}

cat("Best value of k for k-means is",best_k_val," with similarity score as",best_score)
cat("Time taken for the simulation:",Sys.time()-t,"s")
