#' 2D grid graph adjacency matrix

#' This function returns an adjacent matrix of 2D grid graph
#' @param n (m) is an integer indicating grid size
#' @keywords graph grid
#' @export
#' @examples
#' library(igraph)
#' out <- my.grid.adjacent(3,4)
#' image(out)
#' plot(graph.adjacency(out))

my.grid.adjacent <- function(n,m=n){
	N <- my.grid.ad.pre(n)
	M <- my.grid.ad.pre(m)
	gr.dec.prod(N,M)
}

#' @export
my.grid.ad.pre <- function(n){
	ret <- diag(rep(1,n-1))
	ret <- rbind(rep(0,n-1),ret)
	ret <- cbind(ret,rep(0,n))
	ret + t(ret)
}
#' @export
gr.dec.prod <- function(m1,m2){
	n1 <- length(m1[,1])
	n2 <- length(m2[,1])
	n <- n1*n2
	N <- as.matrix(expand.grid(1:n1,1:n2))
	M <- matrix(FALSE,n,n)
	for(i in 1:n){
		for(j in 1:n){
			if((N[i,1] == N[j,1]) & m2[N[i,2],N[j,2]]){
				M[i,j] <- TRUE
			}else if((N[i,2] == N[j,2]) & m1[N[i,1],N[j,1]]){
				M[i,j] <- TRUE
			}
		}
	}
	M
}
