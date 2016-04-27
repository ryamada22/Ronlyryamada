#' Generate Sphere field
#'
#' Sphere field generator using harmonc fx and transformation matrix
#' @param n depth of sphere harmonic functions
#' @param k  simple sphere radius
#' @param n.mesh  fineness of mesh
#' @param scale.shift scale of shift
#' @param scale.rotation1 scale of rotation 1
#' @param scale.rotation2 scale of rotation 2
#' @param scale.shear scale of shear
#' @export
#' @examples
#' library(rgl)
#' library(RFOC)
#' out <- my.sphere.field()
#' # ”¼Œa‚ðF‚Å
#' col1 <- col2 <- rep(0,length(out$v[,1]))
#' col1[which(out$v[,1]>=0)] <- out$v[,1][which(out$v[,1]>=0)]
#' col2[which(out$v[,1]<0)] <- -out$v[,1][which(out$v[,1]<0)]
#' col1 <- col1/max(col1)
#' col2 <- col2/max(col2)
#' plot3d(out$mesh$x,col=rgb(1-col1,1-col2,0),size=10)

my.sphere.field <- function(n=5,k=runif(1),n.mesh=32,scale.shift=0.1,scale.rotation1=0.1,scale.rotation2=0.1,scale.shear=0.1){
	A. <- matrix(runif(n^2),n,n)
	A.[1,1] <- k
	B <- matrix(rnorm(n^2),n,n)
	#B <- matrix(0,n,n)
	xxx <- my.spherical.harm.mesh(A=A.,B=B,n=n.mesh)

	M <- diag(rep(1,4))
	M[1:3,4] <- runif(3)*scale.shift
	theta1 <- runif(1)*scale.rotation1
	theta2 <- runif(1)*scale.rotation2
	Mm <- diag(rep(1,3))
	Mm[1:2,1:2] <- my.2d.rot(theta1) %*% Mm[1:2,1:2]
	Mm[2:3,2:3] <- my.2d.rot(theta2) %*% Mm[2:3,2:3]
	M[1:3,1:3] <- Mm
	M[4,1:3] <- rnorm(3)*scale.shear
	xxxx <- cbind(xxx$v,rep(1,length(xxx$v[,1])))
	rot.xxxx <- t(M %*% t(xxxx))
	rot.xxxx. <- rot.xxxx[,1:3]/rot.xxxx[,4]

	sp.mesh <- my_sphere_tri_mesh(n.mesh)

	return(list(v = rot.xxxx., mesh=sp.mesh))
}
