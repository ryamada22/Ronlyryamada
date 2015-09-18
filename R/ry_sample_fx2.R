#' A Sample Function2
#'
#' This is a toy function
#' @param n is integer for number of randam values.
#' @keywords toy
#' @export
#' @examples
#' n <- 100
#' m <- 3
#' ry_sample_fx2(n,m)

ry_sample_fx2 <- function(n,m){
	plot(rnorm(n,mean=m))
}