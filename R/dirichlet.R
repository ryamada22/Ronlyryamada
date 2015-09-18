#' @export
#' @importFrom MCMCpack rdirichlet

myfxDirichlet <- function(n,m){
	MCMCpack::rdirichlet(n,rep(1,m))
}