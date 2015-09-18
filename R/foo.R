#' @export
#' @useDynLib Rryamada

foo <- function(x) {
    stopifnot(is.numeric(x))

    out <- .Fortran("foo", x = as.double(x), n = length(x), PACKAGE = "Rryamada")
    return(out$x)
}