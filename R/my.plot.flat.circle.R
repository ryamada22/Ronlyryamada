#' place on sphere with radius R
#'
#' Fit to sphere
#' @param x, a matrix; num. of column = dimension
#' @param R is a real value of radius
#' @export
#' @examples
#' x <- matrix(rnorm(100*3),ncol=3)
#' x <- my.fit.sphere(x)

my.fit.sphere <- function(x,R=1){
	r <- sqrt(apply(x^2,1,sum))
	x <- x/r*R
	return(x)
}


#' Plot sphere in a circular disc
#'
#' plot sphere in a curcular disc
#' @param x 2-column matrix with angles
#' @export
#' @examples
#' require(RFOC)
#' x <- matrix(rnorm(100*3),ncol=3)
#' x <- my.fit.sphere(x)
#' x.sp <- TOSPHERE(x[,1],x[,2],x[,3])
#' plot(x.sp[[1]],x.sp[[2]])
#' x.sp.2 <- list(x.sp[[1]]/360*2*pi,x.sp[[2]]/360*2*pi)
#' my.plot.flat.circle(x.sp.2)

my.plot.flat.circle <- function(x){
	plot(x[[2]] * cos(x[[1]]),x[[2]] * sin(x[[1]]))

}


