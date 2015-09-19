#' Weight border points for their local points
#'
#' This function weights points with fraction of points
#' closest to them among all border points
#' @param x is a 0/1 matrix.
#' @param r1 is a positive real value of width to measure circumference.
#' @param r2 is a positive real value to construct border points
#' @keywords shape circumference
#' @export
#' @examples
#' library(igraph)
#' n <- 2^5
#' x <- y <- 1:n
#' xy <- as.matrix(expand.grid(x,y))
#' xy.val <- rep(0,length(xy[,1]))
#' n.circle <- 30
#' ctr <- matrix(sample(20:40,replace=TRUE,n.circle*2),ncol=2)
#' r <- sample(1:20,replace=TRUE,n.circle)
#' 
#' for(i in 1:n.circle){
#' 	tmp.x <- xy[,1] - ctr[i,1]
#' 	tmp.y <- xy[,2] - ctr[i,2]
#' 	s <- which(tmp.x^2+tmp.y^2 < r[i]^2)
#' 	xy.val[s] <- 1
#' }
#' xy.mat <- matrix(xy.val,ncol=n)
#' yyy <- my.dcirc.measure(xy.mat,10,1)
#' bdr <- yyy$border
#' bdr <- rbind(bdr,bdr[2,])
#' diff.bdr <- apply(bdr,2,diff)
#' #diff.bdr <- diff.bdr[-length(diff.bdr[,1]),]
#' diff.bdr.1 <- diff.bdr[-length(diff.bdr[,1]),]
#' diff.bdr.2 <- diff.bdr[-1,]
#' diff.bdr.1.cplx <- diff.bdr.1[,1] + diff.bdr.1[,2]*1i
#' diff.bdr.2.cplx <- diff.bdr.2[,1] + diff.bdr.2[,2]*1i
#' direc <- diff.bdr.2.cplx/diff.bdr.1.cplx
#' arg.direc <- Arg(direc)
#' cumsum.arg <- cumsum(arg.direc)
#' plot(cumsum.arg,type="l")
#' plot(yyy[[1]])
#' plot(yyy[[3]])
#' col <- yyy[[1]]
#' col <- (max(col)-col)/(max(col)-min(col))
#' col <- 0.2 + 0.6*col
#' plot(yyy[[3]],pch=20,col=gray(col))
#' loc <- cumsum(yyy[[1]])
#' plot(loc)
#' plot(loc,cumsum.arg)

my.dcirc.measure <- function(x,r1,r2=1){
	out1 <- my.circ.measure(x,r1)
	bdr <- my.border(x,r2)$border
	bdr. <- bdr[-length(bdr[,1]),]
	bnd <- t(which(out1$band!=0,arr.ind=TRUE))
	len <- matrix(0,length(bnd[1,]),length(bdr.[,1]))
	for(i in 1:length(len[1,])){
		tmp <- (bnd - bdr.[i,])
		tmp2 <- apply(tmp^2,2,sum)
		len[,i] <- tmp2
	}
	min.len <- apply(len,1,min)
	len2 <- len-min.len
	len2 <- len2==0
	len2 <- len2/apply(len2,1,sum)
	len3 <- apply(len2,2,sum)
	return(list(local.length = len3/sum(len3),circ=out1,border=bdr))
}


