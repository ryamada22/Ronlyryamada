#' 2D area's border
#'
#' These function is for 2D 0/1 image matrix and returns its border
#' @param x is a 0/1 matrix.
#' @param r is a positive real value of width to measure circumference.
#' @keywords shape circumference
#' @export
#' @examples
#' library(igraph)
#' n <- 2^5
#' x <- y <- 1:n
#' xy <- as.matrix(expand.grid(x,y))
#' xy.val <- rep(0,length(xy[,1]))
#' n.circle <- 30
#' ctr <- matrix(sample(50:200,replace=TRUE,n.circle*2),ncol=2)
#' r <- sample(1:40,replace=TRUE,n.circle)
#' 
#' for(i in 1:n.circle){
#' 	tmp.x <- xy[,1] - ctr[i,1]
#' 	tmp.y <- xy[,2] - ctr[i,2]
#' 	s <- which(tmp.x^2+tmp.y^2 < r[i]^2)
#' 	xy.val[s] <- 1
#' }
#' xy.mat <- matrix(xy.val,ncol=n)
#' out <- my.border(xy.mat)
#' plot(out$border,type="l")

#library(igraph)

my.border <- function(x,r=2){
	out1 <- my.circ.measure(x,r)
	border.xy <- which(out1$band!=0,arr.ind=TRUE)
	D <- as.matrix(dist(border.xy,method="manhattan"))
	D. <- D==1
	g <- graph.adjacency(D.,mode="undirected")
	mst.g <- minimum.spanning.tree(g)
	D.mst <- get.adjacency(mst.g)
	ends.mst <- farthest_vertices(mst.g)
	end.to.end <- get.shortest.paths(mst.g,ends.mst$vertices[1],ends.mst$vertices[2])
	vs <- as.numeric(end.to.end[[1]][[1]]$name)
	vs <- c(vs,vs[1])
	return(list(v=vs,border=x[vs,]))
}
