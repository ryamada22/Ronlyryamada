#' Schwarz-Christoffel-like curve
#'
#' This function gives Schwarz-Christoffel-like curve.
#' @param p is a positive real number.
#' @param n is a positive integer; No. steps.
#' @param S0 is a complex number with 0 imaginary.
#' @param p2 is a positive real number >= p.
#' @keywords Schwarz-Crhistoffel
#' @export
#' @examples
#' out <- my_schwarzChristoffel(0.2,30)
#' plot(out)

my_schwarzChristoffel <- function(p,n,S0=0+0*1i,p2 = NULL){
	if(is.null(p2)){
		p2 <- p
	}
	rw <- sample(c(-1,1),n,replace=TRUE) * runif(n,min=p,max=p2) + 1/2 + 0 * 1i
	my.S <- function(z,a){
		z <- z + 0*1i
		(z+2*sqrt((1-a)/a))^(1-a)*(z-2*sqrt(a/(1-a)))^a
	}
	S0 <- 0 + 0 * 1i
	Ss <- rep(0,n)

	for(i in 1:n){
		if(i == 1){
			Ss[i] <- my.S(S0,rw[i])
		}else{
			Ss[i] <- S0
			for(j in i:1){
				Ss[i] <- my.S(Ss[i],rw[j])
			}
		}
	}
	Ss
}