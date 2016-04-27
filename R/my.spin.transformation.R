#' Spin transformation
#'
#' Spin transformation
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
#' library(rgl)
#' library(RFOC)
#' library(onion)
#' out <- my.sphere.field()
#' # ”¼Œa‚ðF‚Å
#' col1 <- col2 <- rep(0,length(out$v[,1]))
#' col1[which(out$v[,1]>=0)] <- out$v[,1][which(out$v[,1]>=0)]
#' col2[which(out$v[,1]<0)] <- -out$v[,1][which(out$v[,1]<0)]
#' col1 <- col1/max(col1)
#' col2 <- col2/max(col2)
#' plot3d(out$mesh$x,col=rgb(1-col1,1-col2,0),size=10)

my.spin.transformation <- function(n=5,k=runif(1),n.mesh=32,n.step = 10,scale.shift=0.1,scale.rotation1=0.1,scale.rotation2=0.1,scale.shear=0.1,file.name = "hoge"){
	out1 <- my.sphere.field(n=n,k=k,n.mesh=n.mesh,scale.shift=scale.shift,scale.rotation1=scale.rotation1,scale.rotation2=scale.rotation2,scale.shear=scale.shear)
	out2 <- my.sphere.field(n=n,k=k,n.mesh=n.mesh,scale.shift=scale.shift,scale.rotation1=scale.rotation1,scale.rotation2=scale.rotation2,scale.shear=scale.shear)
	mesh <- out1$mesh
	xyz <- mesh$x
	xyz.h <- xyz[,1] * Hi + xyz[,2] * Hj + xyz[,3] * Hk
	rot.h <- out2$v[,1] + out1$v[,1] * Hi + out1$v[,2] * Hj + out1$v[,3] * Hk
	new.xyz.h <- Conj(rot.h) * xyz.h * rot.h
	new.xyz <- cbind(i(new.xyz.h),j(new.xyz.h),k(new.xyz.h))
	return(list(new.xyz=new.xyz,mesh=mesh))
}
