#' Arbitrary moment
#'
#' This function returns order-th weighted sample moment.
#' @param x is a numeric vector.
#' @param p is a numeric vector of weight of x elements
#' @param order is order of moment.
#' @param center is logical; when TRUE, moment is around mean and otherwise zero.
#' @keywords moment
#' @export
#' @examples
#' n <- 10
#' x <- runif(n)
#' p <- runif(n)
#' p <- p/sum(p)
#' out <- ry_momentProb(x,p,oder=2,center=TRUE)

ry_momentProb<-function (x, p,order = 1, center = FALSE) 
{
 if(center)
  x <- x - ry_momentProb(x,p,order=1,center=FALSE)
 sum(x^order*p)
}

