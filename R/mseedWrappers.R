# wrapper function to invoke test_hello_mseed
test_hello_mseed <- function() {
  result <- .Call("test_hello_mseed")
  return(result)
}