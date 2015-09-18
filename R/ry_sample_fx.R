#' A Sample Function
#'
#' This is a toy function
#' @param n is integer for number of randam values.
#' @keywords toy
#' @export
#' @examples
#' n <- 100
#' ry_sample_fx(n)

ry_sample_fx <- function(n){
	plot(rnorm(n))
}