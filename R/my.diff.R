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
#' my.diff.serial(pcaout.slice[[5]][[2]])

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

#' @export

my.arr.address <- function(addr,dm){
	tmp <- cumprod(c(1,dm))
	(addr-1) %*% tmp[-length(tmp)]+1
}
	
#' @export

my.arr.address.inv <- function(v,dm){
	tmp <- cumprod(c(1,dm))
	v. <- v-1
	ret <- matrix(0,length(v),length(dm))
	for(i in 1:length(dm)){
		ret[,i] <- (v. %% tmp[i+1]) / tmp[i]
		v. <- v. - (ret[,i] * tmp[i])
	}
	ret + 1
}

#' @export

my.shift <- function(v,u,dm){
	addr <- 1:prod(dm)
	arr.addr <- my.arr.address.inv(addr,dm)
	shifted.addr <- t(t(arr.addr) + u)
	addr.within <- which(apply((shifted.addr>0) & (t(t(shifted.addr)- dm)<=0),1,prod)==1)
	
	shifted.addr <- my.arr.address(shifted.addr[addr.within,],dm)
	shifted.v <- rep(0,prod(dm))
	shifted.v[shifted.addr] <- v[addr.within]
	return(shifted.v)	
}

#' @export

my.shift.patterns <- function(x){
	x.floor <- floor(x)
	d <- length(x)
	tmp <- list()
	for(i in 1:d){
		tmp[[i]] <- c(0,1)
	}
	tmp2 <- as.matrix(expand.grid(tmp))
	ret.x <- t(t(tmp2) + x.floor)
	
	p <- x-x.floor
	p. <- 1-p
	tmp <- list()
	for(i in 1:d){
		tmp[[i]] <- c(p[i],p.[i])
	}
	tmp2 <- as.matrix(expand.grid(tmp))
	ret.p <- apply(tmp2,1,prod)
	
	return(list(x = ret.x,p=ret.p))
}

#' @export

my.diff.shifted <- function(v1,v2,shift){
	shift.patterns <- my.shift.patterns(shift)
	ret <- v1
	for(i in 1:length(shift.patterns[[2]])){
		tmp <- my.shift(c(v2),shift.patterns[[1]][i,],dim(v2))
		shifted.v2 <- array(tmp,dim(v2))
		ret <- ret - shift.patterns[[2]][i] * shifted.v2
	}
	return(ret)
}

#' @export

my.diff.center.shift <- function(v1,v2){
	shift <- my.whole.center(v1) - my.whole.center(v2)
	my.diff.shifted(v1,v2,shift)
}