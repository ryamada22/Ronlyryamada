#' e-flat and m-flat parameterizations for {0,1}^d
#'
#' e-flat and m-flat
#' @export
#' @examples
#' library(MCMCpack)
#' d <- 3
#' n <- 2^d
#' freq <- c(rdirichlet(1,rep(1,n)))
#' freq
#' thetapsi <- my.freq2thetapsi(freq)
#' psi <- thetapsi[1]
#' theta <- thetapsi[-1]
#' thetapsi
#' eta <- my.freq2eta(freq)
#' eta
#' my.eta2freq(eta)
#' my.eta2thetapsi(eta)
#' my.theta2freq(theta)
#' my.theta2eta(theta)
#' 
#' freq <- rep(1/8,8)
#' freq
#' thetapsi <- my.freq2thetapsi(freq)
#' psi <- thetapsi[1]
#' theta <- thetapsi[-1]
#' thetapsi
#' eta <- my.freq2eta(freq)
#' eta
#' my.eta2freq(eta)
#' my.eta2thetapsi(eta)
#' my.theta2freq(theta)
#' my.theta2eta(theta)
#' 
#' dd <- 0
#' freq <- c(1/2-dd*3,dd,dd,dd,dd,dd,dd,1/2-dd*3)
#' freq
#' thetapsi <- my.freq2thetapsi(freq)
#' psi <- thetapsi[1]
#' theta <- thetapsi[-1]
#' thetapsi
#' eta <- my.freq2eta(freq)
#' eta
#' my.eta2freq(eta)
#' my.eta2thetapsi(eta)
#' my.theta2freq(theta)
#' my.theta2eta(theta)

my.haplotypes <- function(d){
	n <- 2^d
	X <- matrix(0,d,n)
	cnt <- 1
	for(i in 1:d){
		tmp <- combn(d,i)
		for(j in 1:length(tmp[1,])){
			X[tmp[,j],cnt+j] <- 1
		}
		cnt <- cnt + length(tmp[1,])
	}
	return(X)
}
#' @export

my.freq2eta.mat <- function(d){
	X <- my.haplotypes(d)
	n <- 2^d
	M <- matrix(0,n,n)

	for(i in 1:n){
		M[,i] <- apply(X >= X[,i],2,prod)
	}
	return(t(M[2:n,2:n]))
}
#' @export

my.theta2freq.mat <- function(d){
	K <- my.freq2eta.mat(d)
	n <- 2^d
	M <- matrix(0,n,n)
	M[2:n,2:n] <- K
	M[,1] <- -1
	return(M)
}
#' @export

my.freq2theta.mat <- function(d){
	M <- my.theta2freq.mat(d)
	return(solve(M))
}
#' @export

my.theta2psi <- function(theta){
	len <- length(theta)
	d <- log(len+1,2)
	M <- my.theta2freq.mat(d)
	#tmp <- M %*% c(0,theta)
	tmp <- my.Inf.matrix.mult(M,c(0,theta))
	delta <- sum(exp(tmp))
	psi <- log(delta)
	return(list(psi=psi,M=M))
}
#' @export

my.theta2freq <- function(theta){
	out <- my.theta2psi(theta)
	psi.theta <- c(out$psi,theta)
	M <- out$M
	return(exp(my.Inf.matrix.mult(M,psi.theta)))
}
#' @export

my.freq2thetapsi <- function(freq){
	d <- log(length(freq),2)
	M <- my.freq2theta.mat(d)
	#return(M %*% log(freq))
	return(my.Inf.matrix.mult(M,log(freq)))
}
#' @export

my.freq2eta <- function(freq){
	d <- log(length(freq),2)
	K <- my.freq2eta.mat(d)
	return(K %*% freq[-1])
}
#' @export

my.eta2freq <- function(eta){
	d <- log(length(eta)+1,2)
	K.inv <- solve(my.freq2eta.mat(d))
	tmp <- K.inv %*% eta
	return(c(1-sum(tmp),tmp))
}
#' @export

my.theta2eta <- function(theta){
	freq <- my.theta2freq(theta)
	return(my.freq2eta(freq))
}
#' @export

my.eta2thetapsi <- function(eta){
	freq <- my.eta2freq(eta)
	return(my.freq2thetapsi(freq))
}
#' @export

my.matrix.Inf.separator <- function(M){
	if(!is.matrix(M)){
		M <- matrix(M,ncol=1)
	}
	M.inf <- M.noninf <- matrix(0,length(M[,1]),length(M[1,]))
	infs <- which(is.infinite(M))
	noninfs <- which(!is.infinite(M))
	M.noninf[noninfs] <- M[noninfs]
	M.inf[infs] <- sign(M[infs])
	return(list(M.noninf=M.noninf,M.inf=M.inf))
}
#' @export

my.Inf.matrix.mult <- function(M1,M2){
	out1 <- my.matrix.Inf.separator(M1)
	out2 <- my.matrix.Inf.separator(M2)
	
	tmp.n.n <- out1$M.noninf %*% out2$M.noninf
	tmp.i.i <- out1$M.inf %*% out2$M.inf
	tmp.n.i <- out1$M.noninf %*% out2$M.inf
	tmp.i.n <- out1$M.inf %*% out2$M.noninf
	
	tmp.any.i <- sign(tmp.i.i + tmp.n.i + tmp.i.n)
	non.zero <- which(tmp.any.i!=0)
	tmp.any.i[non.zero] <- tmp.any.i[non.zero] * Inf
	tmp <- tmp.n.n + tmp.any.i
	return(tmp)
}
