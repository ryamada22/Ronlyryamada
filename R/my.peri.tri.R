#' 3D object surface triangles
#'
#' These function returns surface trinangles of 3D object
#' @param v is a 3d-array..
#' @param h is a real value to draw contour line.
#' @keywords shape triangulation mesh
#' @export
#' @examples
#' library(misc3d)
#' k <- 4
#' out <- my.3d.obj(3,k)
#' obj <- out$xyz[which(out$xyz.val==1),]
#' library(rgl)
#' plot3d(obj)
#' 
#' xyz.arr <- array(out$xyz.val,c(2^k,2^k,2^k))
#' tris <- my.peri.tri(xyz.arr)
#' drawScene.rgl(tris)

my.peri.tri <- function(v,h=max(v)){
	#library(misc3d)
	con <- computeContour3d(v, h, 1)
	tris <- makeTriangles(con)

	V1 <- tris$v1
	V2 <- tris$v2
	V3 <- tris$v3

	L12 <- apply((V1-V2)^2,1,sum)
	L13 <- apply((V1-V3)^2,1,sum)
	L123 <- L12+L13
	zeros <- (L123==0)

	tris$v1 <- V1[!zeros,]
	tris$v2 <- V2[!zeros,]
	tris$v3 <- V3[!zeros,]
	tris	
}

#' @export

my.3d.obj <- function(n,k=5,a=6,b=4){
	x <- 1:(2^k)
	xyz <- as.matrix(expand.grid(x,x,x))
	xyz.val <- rep(0,length(xyz[,1]))
	rs <- runif(n) * (max(x)/a)
	ctrs <- matrix(runif(n*3),ncol=3) * max(x)/b + max(x)*(b-1)/(b*2)
	for(i in 1:n){
		tmp <- which(apply((t(xyz) - ctrs[i,])^2,2,sum) < rs[i]^2)
		xyz.val[tmp] <- 1
	}
	return(list(xyz=xyz,xyz.val=xyz.val))
}

