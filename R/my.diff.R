#' Difference evaluation of scalers/vectors/matrices/arrays/others

#' The functions are to return difference between two consecutive data objects
#' Difference of the 1st axis should be taken so that simple vector handling is possible.
#' @export
#' @examples
#' library(igraph)
#' library(misc3d)
#' library(onion)
#' library(rgl)
#' df <- 4
#' n <- 20
#' data. <- array(0,rep(n,df))
#' addr <- t(which(data.>=0,arr.ind=TRUE))
#' n.sp <- 6
#' for(i in 1:n.sp){
#' 	tmp <- runif(df) * n
#' 	tmp2 <- runif(1) * n/5+2
#' 	data.[apply((addr - tmp)^2,2,sum) < tmp2] <- 1
#' }
#' clu <- my.labeling(data.)
#' hist(clu[[3]],ylim=c(0,100))
#' cluster.n <- clu[[1]][[3]]
#' vol <- vol.slice <- list()
#' ctr <- ctr.slice <- list()
#' #tri.mesh.slice <- list()
#' pcaout <- pcaout.slice <- list()
#' dm <- dim(data.)
#' for(i in 1:cluster.n){
#' 	tmp.arr <- array(as.numeric(clu[[3]]==i),dm)
#' 	vol[[i]] <- my.whole.vol(tmp.arr)
#' 	vol.slice[[i]] <- my.slice.vol(tmp.arr,df)
#' 	ctr[[i]] <- my.whole.center(tmp.arr)
#' 	ctr.slice[[i]] <- my.slice.center(tmp.arr,df)
#' 	pcaout[[i]] <- my.whole.pca(tmp.arr)
#' 	pcaout.slice[[i]] <- my.slice.pca(tmp.arr,df)
#' }
#' my.diff.serial(ctr.slice[[5]])
#' my.diff.serial(pcaout.slice[[5]])

my.diff <- function(v){
	if(is.vector(v)){
		v <- matrix(v,ncol=1)
	}
	dm <- dim(v)
	n <- prod(dm[-length(dm)])
	return(array(c(v[(n+1):length(v)]) - c(v[1:(length(v)-n)]),c(dm[-length(dm)],dm[length(dm)]-1)))
}

#' @export

my.diff.serial <- function(v,k=1){
	ret <- list()
	ret[[1]] <- my.diff(v)
	if(k>1){
		for(i in 2:k){
			ret[[i]] <- my.diff(ret[[i-1]])
		}
	}
	return(ret)
}
