#' Arbitrary dimensional voxel data labeling

#' This function labels continuous segments.
#' @param data is a matrix with ncol-dimensions
#' @export
#' @examples
#' library(igraph)
#' df <- 4
#' n <- 20
#' data <- array(0,rep(n,df))
#' addr <- t(which(data>=0,arr.ind=TRUE))
#' n.sp <- 20
#' for(i in 1:n.sp){
#' 	tmp <- runif(df) * n
#' 	tmp2 <- runif(1) * n/5
#' 	data[apply((addr - tmp)^2,2,sum) < tmp2] <- 1
#' }
#' clu <- my.labeling(data)
#' hist(clu[[3]],ylim=c(0,100))
#' dm <- c(2^8,2^8,2^3,2^4)
#' ar <- array(1:prod(dm),dm)
#' d <- sample(1:length(dm),1)
#' j <- sample(1:dm[d],1)
#' out1 <- my.slice(c(ar),dm,d,j)
#' out2 <- my.slice.2(c(ar),dm,d,j)
#' plot(out2[[1]])

my.labeling <- function(data){
	x <- which(data==1,arr.ind=TRUE)
	d <- as.matrix(dist(x,method="manhattan"))
	d <- d==1
	g <- graph.adjacency(d)
	clu <- components(g)
	addr <- which(data>=0,arr.ind=TRUE)
	val <- rep(0,length(addr[,1]))
	dm <- dim(data)
	dm.prod <- c(1,cumprod(dm))
	dm.prod <- dm.prod[-length(dm.prod)]
	for(i in 1:clu[[3]]){
		tmp <- matrix(x[clu[[1]]==i,],ncol=ncol(x))
		tmp2 <- (tmp-1) %*% dm.prod + 1
		val[tmp2] <- i
	}
	return(list(clu=clu,dm=dm,val=val))
	#groups(clu)
}

#' @export
my.slice <- function(v,dm,d,j){
	tmp <- list()
	for(i in 1:length(dm)){
		tmp[[i]] <- 1:dm[i]
	}
	tmp2 <- as.matrix(expand.grid(tmp))
	return(list(v=v[which(tmp2[,d]==j)],dm=dm[-d]))
}

#' @export

my.slice.2 <- function(v,dm,d,j){
	j <- j-1
	n <- length(dm)
	dm. <- c(1,dm)
	n <- length(dm.)
	pre <- dm.[1:d]
	post <- dm.[1:(d+1)]
	pre.prod <- prod(pre)
	post.prod <- prod(post)
	all.prod <- prod(dm)
	n.rep <- all.prod/post.prod
	tmp1 <- (j)*pre.prod + (1:pre.prod)
	tmp2 <- ((0:(n.rep-1))) * post.prod
	addr <- outer(tmp1,tmp2,"+") 
	return(list(v=v[addr],dm=dm[-d]))
}
