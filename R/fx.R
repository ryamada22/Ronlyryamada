#' @export
myfx1 <- function(n){
	runif(n)
}
#' @export
myfx2 <- function(n){
	rnorm(n)
}

myfx3 <- function(){
	sample(1:100,1)
}

#' @export
myfx4 <- function(){
	n <- myfx3()
	rnorm(n)
}