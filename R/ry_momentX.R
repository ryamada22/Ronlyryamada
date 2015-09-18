#' Arbitrary moment
#'
#' This function returns order-th sample moment.
#' @param x is a numeric vector.
#' @param order is order of moment.
#' @param center is logical; when TRUE, moment is around mean otherwise zero
#' @keywords complex function, conformal transformation
#' @export
#' @examples
#' n <- 10
#' x <- runif(n)
#' out2 <- ry_momentX(x,order=2,center=TRUE)


ry_momentX<-function (x, order = 1, center = FALSE) 
{
 ry_momentProb(x,p=rep(1,length(x))/length(x),order=order,center=center) 
}