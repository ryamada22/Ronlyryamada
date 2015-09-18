#' Mobius transformation
#'
#' This function returns a moebius transformation
#' function that maps three complex points, zs, to
#' three complex points, ws.
#' @param zs is a complex vector with length three.
#' @param ws is a complex vector with length three.
#' @keywords complex function, conformal transformation
#' @export
#' @examples
#' zs <- runif(3) + 1i * runif(3)
#' ws <- runif(3) + 1i * runif(3)
#' mob.fx <- ry_moebius(zs,ws)
#' print(mob.fx(zs))
#' print(ws)


ry_moebius <- function(zs,ws){
	z1 <- zs[1];z2 <- zs[2];z3 <- zs[3];
	w1 <- ws[1];w2 <- ws[2];w3 <- ws[3];
	ret <- function(z){
		(w1*(w3-w2)*(z3-z1)*(z-z2)-w2*(w3-w1)*(z3-z2)*(z-z1))/((w3-w2)*(z3-z1)*(z-z2)-(w3-w1)*(z3-z2)*(z-z1))
	}
	ret
}