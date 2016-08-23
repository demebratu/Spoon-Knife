uthor: Aaron Gabriel
# 8/3/16
# K Harmonic Means R Code

# clear workspace
rm(list=ls())

# initialize data
x <- matrix( nrow=7,ncol=2, c(1.0,1.5,3.0,5.0,3.5,4.5,3.5,1.0,2.0,4.0,7.0,5.0,5.0,4.5)  )

# initialize centroids
centroids <- matrix(nrow=2,ncol=2,c(1.0,5.0,1.0,7.0))
 
khmeans <- function(x, centroids, maxIterations) {
  
  numClusters <- nrow(centroids)
  
  numDataPoints <- nrow(x)
  
  d <- matrix(0, numClusters, numDataPoints)
  
  dharmonic <- d
  
  labels <- rep(0, numDataPoints)
  
  for (run in 1:maxIterations) {
  
    for (i in 1:numClusters) {
	
      for (j in 1:numDataPoints) {
        d[i,j] <- sqrt(sum((centroids[i,] - x[j,])^2))  # d[i,j] - distance of centroid i to data point j
      }
	  
    }
	
    labels <- apply(d, 2, which.min) # update clusters
    d <- d + 1e-6                    # very easy to get NaN with this implementation
        	
    for (i in 1:numClusters) {
	
      for (j in 1:numDataPoints) {
	  
        dharmonic[i,j] <- (d[i,j]^3)*((sum(1/d[,j]^2))^2)      # I will use this form of harmonic mean for the
		                                                       # weights to use in centroid update, so it will look
      }                                                        # like a dot product

	}  
       
    for (i in 1:numClusters) {      
	
      centroids[i,] <- colSums(dharmonic[i,labels == i]*x[labels == i,]) / sum(dharmonic[i,labels == i])  # to update centroid i, get my x's in those clusters and dot product each dimension with
                                                                                                          # the row formed by taking ith row of my harmonics matrix (only columns corresponding to that cluster)
	  }                                                                                                   # divide by sum of those harmonics corresponding to that cluster
                                                                                                          
  }
  
  return(list(labels, centroids))  
}
