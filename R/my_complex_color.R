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
#' out <- my_complex_color(3+1i*2)

my_complex_color <- function(z,int0=0.6,sat0=0.3,int1=1,sat1=1){
my.hsv <- function(z,int0=0.6,sat0=0.3,int1=1,sat1=1){
arg <- Arg(z)
s <- which(arg<0)
arg[s] <- arg[s]+2*pi
r <- Mod(z)
s <- which(r>1)
r[s] <- log(r[s])
r. <- 4*(r%%1)
k <- floor(r.)
r. <- r.-k
inten <- sat <- rep(0,length(r))
s <- which(k==0)
inten[s] <- int1
sat[s] <- sat1-(sat1-sat0)*r.[s]
s <- which(k==1)
inten[s] <- int1-(int1-int0)*r.[s]
sat[s] <- sat0
s <- which(k==2)
inten[s] <- int0
sat[s] <- sat1-(sat1-sat0)*(1-r.[s])
s <- which(k==3)
inten[s] <- int1-(int1-int0)*(1-r.[s])
sat[s] <- sat1

return(cbind(arg,inten,sat))
}
my.hsv2rgb <- function(h,s,v){
hi <- floor(h/(2*pi)*6)
hi[which(hi==6)] <- 0
f <- (h/(2*pi)*6) %%1
p <- v*(1-s)
q <- v *(1-f*s)
t <- v *(1-(1-f)*s)
r <- g <- b <- rep(0,length(h))
s <- which(hi==0)
r[s] <- v[s];g[s] <- t[s]; b[s] = p[s];
s <- which(hi==1)
r[s] <- q[s];g[s] <- v[s]; b[s] = p[s];
s <- which(hi==2)
r[s] <- p[s];g[s] <- v[s]; b[s] = t[s];
s <- which(hi==3)

r[s] <- p[s];g[s] <- q[s]; b[s] = v[s];
s <- which(hi==4)
r[s] <- t[s];g[s] <- p[s]; b[s] = v[s];
s <- which(hi==5)
r[s] <- v[s];g[s] <- p[s]; b[s] = q[s];
return(cbind(r,g,b))
}
hsv <- my.hsv(z,int0=int0,sat0=sat0,int1=int1,sat1=sat1)
col <- my.hsv2rgb(hsv[,1],hsv[,3],hsv[,2])
rgb(col[,1],col[,2],col[,3])

}