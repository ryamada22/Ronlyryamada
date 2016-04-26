#' place on sphere with radius R
#'
#' Fit to sphere
#' @param x, a matrix; num. of column = dimension
#' @param R is a real value of radius
#' @export
#' @examples
#' x <- matrix(rnorm(1000*3),ncol=3)
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
#' x <- matrix(rnorm(1000*3),ncol=3)
#' x <- my.fit.sphere(x)
#' x.sp <- TOSPHERE(x[,1],x[,2],x[,3])
#' plot(x.sp[[1]],x.sp[[2]])
#' x.sp.2 <- list(x.sp[[1]]/360*2*pi,x.sp[[2]]/360*2*pi)
#' my.plot.flat.circle(x.sp.2)

my.plot.flat.circle <- function(x){
	plot(x[[2]] * cos(x[[1]]),x[[2]] * sin(x[[1]]))

}

#' Plot two hemi-spheres 
#'
#' plot two hemi-spheres
#' @param x 2-column matrix with angles
#' @export
#' @examples
#' require(RFOC)
#' x <- matrix(rnorm(1000*3),ncol=3)
#' x <- my.fit.sphere(x)
#' x.sp <- TOSPHERE(x[,1],x[,2],x[,3])
#' plot(x.sp[[1]],x.sp[[2]])
#' x.sp.2 <- list(x.sp[[1]]/360*2*pi,x.sp[[2]]/360*2*pi)
#' my.plot.hemispheres(x.sp.2)

my.plot.hemispheres <- function(x){
	n <- which(x[[2]] <= pi/2)
	s <- which(x[[2]] > pi/2)
	
	N <- cbind(x[[2]][n] * cos(x[[1]][n]),x[[2]][n] * sin(x[[1]][n]))
	s.r <- pi-x[[2]][s]
	S <- cbind(s.r * cos(x[[1]][s]), s.r * sin(x[[1]][s]))
	S[,2] <- S[,2] - pi
	NS <- rbind(N,S)
	plot(NS,xlim=range(NS),ylim=range(NS))

}
