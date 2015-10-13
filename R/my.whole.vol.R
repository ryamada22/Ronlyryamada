#' Featuring 3D+1T array data

#' The functions are to quantitate various features of 3D+1T data
#' @export
#' @examples
#' library(igraph)
#' df <- 4
#' n <- 20
#' data. <- array(0,rep(n,df))
#' addr <- t(which(data.>=0,arr.ind=TRUE))
#' n.sp <- 20
#' for(i in 1:n.sp){
#' 	tmp <- runif(df) * n
#' 	tmp2 <- runif(1) * n/5
#' 	data.[apply((addr - tmp)^2,2,sum) < tmp2] <- 1
#' }
#' clu <- my.labeling(data.)
#' hist(clu[[3]],ylim=c(0,100))
#' cluster.n <- clu[[1]][[3]]
#' vol <- vol.slice <- list()
#' ctr <- ctr.slice <- list()
#' tri.mesh.slice <- list()
#' dm <- dim(data.)
#' for(i in 1:cluster.n){
#' 	tmp.arr <- array(as.numeric(clu[[3]]==i),dm)
#' 	vol[[i]] <- my.whole.vol(tmp.arr)
#' 	vol.slice[[i]] <- my.slice.vol(tmp.arr,df)
#' 	ctr[[i]] <- my.whole.center(tmp.arr)
#' 	ctr.slice[[i]] <- my.slice.center(tmp.arr,df)
#' 	tri.mesh.slice[[i]] <- my.slice.tri.mesh(tmp.arr.df)
#' }

my.whole.vol <- function(v){
	sum(v)
}

#' @export
my.slice.vol <- function(v,d){
	dm <-dim(v)
	L <- dm[d]
	ret <- rep(NA,L)
	for(i in 1:L){
		ret[i] <- my.whole.vol(my.slice.2(c(v),dm,d,i)[[1]])
	}
	return(ret)
}

#' @export
my.whole.center <- function(v){
	if(!is.matrix(v) & !is.array(v)){
		v <- matrix(v,nrow=1)
	}
	tmp <- which(v==1,arr.ind=TRUE)
	return(apply(tmp,2,mean))
}

#' @export
my.slice.center <- function(v,d){
	dm <- dim(v)
	L <- dm[d]
	ret <- matrix(NA,L,length(dm))
	for(i in 1:L){
		tmp.out <- my.slice.2(c(v),dm,d,i)
		tmp <- my.whole.center(array(tmp.out[[1]],tmp.out[[2]]))
		if(d==1){
			tmp <- c(i,tmp)
		}else if(d==length(dm)){
			tmp <- c(tmp,i)
		}else{
			pre <- 1:(d-1)
			post <- d:length(dm)
			tmp2 <- c(tmp[pre],i,tmp[post])
			tmp <- tmp2
		}
		ret[i,] <- tmp
	}
	return(ret)
}

#' @export
my.slice.tri.mesh <- function(v,d){
	dm <- dim(v)
	L <- dm[d]
	ret <- list()
	for(i in 1:L){
		tmp.out <- my.slice.2(c(v),dm,d,i)
		ret[[i]] <- my.peri.tri(array(tmp.out[[1]],tmp.out[[2]]))
	}
	return(ret)
}