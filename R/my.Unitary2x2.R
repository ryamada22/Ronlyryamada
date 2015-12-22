#' 2x2 Unitary Matrix with determinant +/- 1
#'
#' 2x2 Unitary Matrix
#' @param theta, a real value.
#' @param phi1, a real value.
#' @param phi2, a real value.
#' @keywords Unitary matrix
#' @export
#' @examples
#' theta <- pi/6;phi1 <- pi/4;phi2 <- pi/3
#' U <- my.Unitary2x2(theta,phi1,phi2)
#' U

my.Unitary2x2 <- function(theta,phi1,phi2){
	x <- cos(theta) * exp(1i*phi1)
	y <- sin(theta) * exp(1i*phi2)
	ret1 <- matrix(c(x,y,-Conj(y),Conj(x)),byrow=TRUE,2,2)
	ret2 <- matrix(c(x,y,Conj(y),-Conj(x)),byrow=TRUE,2,2)
	return(list(u1 = ret1,u2=ret2))
}

#' @export

my.complexV2 <- function(theta,phi1,phi2){
	ret <- c(cos(theta)*exp(1i*phi1),sin(theta)*exp(1i*phi2))
	return(ret)
}
