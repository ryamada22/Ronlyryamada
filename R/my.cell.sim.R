#' Cell transformation simulation
#'
#' Simulate cell transformation and write series of obj files
#' @param n depth of sphere harmonic functions
#' @param k  simple sphere radius
#' @param n.mesh  fineness of mesh
#' @param n.step number of time steps
#' @param scale.shift scale of shift
#' @param scale.rotation1 scale of rotation 1
#' @param scale.rotation2 scale of rotation 2
#' @param scale.shear scale of shear
#' @param file.name output file name
#' @export
#' @examples
#' library(devtools)
#' install_github("ryamada22/Ronlyryamada")
#' library(Ronlyryamada)
#' library(rgl)
#' library(RFOC)
#' # parameters
#' n <- 5 # depth of sphere harmonic functions
#' k <- 5 # simple sphere radius
#' n.mesh <- 32 # fineness of mesh
#' n.step <- 50 # number of time steps
#' scale.shift <- 0.1 # shift
#' scale.rotation1 <- 0.1 # rotation 1
#' scale.rotation2 <- 0.1 # rotation 2
#' scale.shear <- 0.1 # shear
#' file.name <- "hoge"
#' my.cell.sim()

my.cell.sim <- function(n=5,k=5,n.mesh=32,n.step=50,scale.shift=0.1,scale.rotation1=0.1,scale.rotation2=0.1,scale.shear=0.1,file.name="hoge"){
	A. <- matrix(runif(n^2),n,n)
	A.[1,1] <- k
	B <- matrix(rnorm(n^2),n,n)
	#B <- matrix(0,n,n)
	xxx <- my.spherical.harm.mesh(A=A.,B=B,n=n.mesh)
	plot3d(xxx$v)
	segments3d(xxx$v[c(t(xxx$edge)),])



	M <- diag(rep(1,4))
	M[1:3,4] <- runif(3)*scale.shift
	theta1 <- runif(1)*scale.rotation1
	theta2 <- runif(1)*scale.rotation2
	Mm <- diag(rep(1,3))
	Mm[1:2,1:2] <- my.2d.rot(theta1) %*% Mm[1:2,1:2]
	Mm[2:3,2:3] <- my.2d.rot(theta2) %*% Mm[2:3,2:3]
	M[1:3,1:3] <- Mm
	M[4,1:3] <- rnorm(3)*scale.shear

	for(i in 1:n.step){
		A. <- A. + rnorm(n^2,0,0.05)
		xxx <- my.spherical.harm.mesh(A=A.,n=32)
		xxxx <- cbind(xxx$v,rep(1,length(xxx$v[,1])))
		rot.xxxx <- t(M %*% t(xxxx))
		rot.xxxx. <- rot.xxxx[,1:3]/rot.xxxx[,4]
		plot3d(rot.xxxx.)
		segments3d(rot.xxxx.[c(t(xxx$edge)),])
		file.out <- paste(file.name,i,".obj",sep="")
		my.write.obj(rot.xxxx.,xxx$f,file.out)
	}


}

#' @export
my.2d.rot <- function(theta){
	matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
}
#' @export
my.write.obj <- function(v,tri,filename){
	fileConn <- file(filename,"w")
	for(i in 1:length(v[,1])){
		tmp <- paste("v",v[i,1],v[i,2],v[i,3],"\n",sep=" ")
		cat(tmp,file=fileConn)
	}
	for(i in 1:length(tri[,1])){
		tmp <- paste("f",tri[i,1],tri[i,2],tri[i,3],"\n",sep=" ")
		cat(tmp,file=fileConn)
	}
	close(fileConn)
}