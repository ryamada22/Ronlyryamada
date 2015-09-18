#' All possible 2x3 tables

#' This function returns nx6 matrix; each row has 6 elements of 2x3 table
#' sharing marginal counts with argument m.
#' @param m is a vector with 6 non-negative values, the first 3 values for cases
#' @keywords SNPs exact trend
#' @export
#' @examples
#' m <- c(9,42,49,6,43,51)
#' my.Cochran.Armitage.trend.exact.2x3(m)
#' my.trend.exact.2x3(m)

my.Cochran.Armitage.trend.exact.2x3 <- function(m){
	out1 <- my.all.2x3(m)
	# marginal counts
	sum.col <- c(m[1]+m[4],m[2]+m[5],m[3]+m[6])
	sum.row <- c(sum(m[1:3]),sum(m[4:6]))
	sum.all <- sum(m)

	tmp1 <-sum(lgamma(sum.col+1))+sum(lgamma(sum.row+1))-lgamma(sum.all+1)
	tmp2 <- apply(lgamma(out1+1),1,sum)
	tmp3 <- tmp1 - tmp2
	
	# Cochran-Armitage
	# http://aoki2.si.gunma-u.ac.jp/R/Cochran-Armitage.html 
	# http://d.hatena.ne.jp/ryamada22/20070601 
	r <- out1[,1:3]
	rr <- out1[,1]*2 + out1[,2]
	rr.range <- range(rr)
	rr.seq <- seq(from=rr.range[1],to=rr.range[2],by=1)
	tmp.pr <- rep(0,length(rr.seq))
	for(i in 1:length(tmp.pr)){
		selec <- which(rr==rr.range[1]+i-1)
		tmp.pr[i] <- sum(exp(tmp3[selec]))
	}
	ori.selec <- which(rr.seq==m[1]*2+m[2])
	#ori.pr <- sum(exp(tmp3[ori.selec]))
	ori.pr <- tmp.pr[ori.selec]
	sum(tmp.pr[which(tmp.pr<=ori.pr)])
}

#' @export
my.Cochran.Armitage.trend.exact.2x3.ori <- function(m){
	out1 <- my.all.2x3(m)
	# marginal counts
	sum.col <- c(m[1]+m[4],m[2]+m[5],m[3]+m[6])
	sum.row <- c(sum(m[1:3]),sum(m[4:6]))
	sum.all <- sum(m)

	tmp1 <-sum(lgamma(sum.col+1))+sum(lgamma(sum.row+1))-lgamma(sum.all+1)
	tmp2 <- apply(lgamma(out1+1),1,sum)
	tmp3 <- tmp1 - tmp2
	
	# Cochran-Armitage
	# http://aoki2.si.gunma-u.ac.jp/R/Cochran-Armitage.html 
	# http://d.hatena.ne.jp/ryamada22/20070601 
	r <- out1[,1:3]
	n <- out1[,1:3] + out1[,4:6]
	R <- sum.row[1]
	N <- sum(m)
	S <- sum.row[2]
	up <- (N-1)*(N*(r[,2]+2*r[,3])-R*(n[,2]+2*n[,3]))^2
	dn <- R*S*(N*(n[,2]+4*n[,3])-(n[,2]+2*n[,3])^2)
	
	out3 <- up/dn
	r.ori <- m[1:3]
	n.ori <- m[1:3] + m[4:6]
	up.ori <- (N-1)*(N*(r.ori[2]+2*r.ori[3])-R*(n.ori[2]+2*n.ori[3]))^2
	dn.ori <- R*S*(N*(n.ori[2]+4*n.ori[3])-(n.ori[2]+2*n.ori[3])^2)
	original.chisq <- up.ori/dn.ori
	selected <- which(abs(out3)>=abs(original.chisq))
	sum(exp(tmp3[selected]))
}

#' @export
my.all.2x3 <- function(m){
	sum.col <- c(m[1]+m[4],m[2]+m[5],m[3]+m[6])
	sum.row <- c(sum(m[1:3]),sum(m[4:6]))
	ceil11 <- min(sum.col[1],sum.row[1])
	ceil12 <- min(sum.col[2],sum.row[1])
	v1 <- 0:ceil11
	v2 <- 0:ceil12
	v12 <- as.matrix(expand.grid(v1,v2))
	v123 <- cbind(v12,sum.row[1]-apply(v12,1,sum))
	v123456 <- cbind(v123,sum.col[1]-v123[,1],sum.col[2]-v123[,2],sum.col[3]-v123[,3])
	v123456. <- v123456[which(apply(v123456,1,min)>=0),]
	v123456.
}


#' @export

my.exact.prob.2x3 <- function(m,log=FALSE){
	sum.col <- c(m[1]+m[4],m[2]+m[5],m[3]+m[6])
	sum.row <- c(sum(m[1:3]),sum(m[4:6]))
	sum.all <- sum(m)
	ret <- sum(lgamma(sum.col+1))+sum(lgamma(sum.row+1))-sum(lgamma(sum.all+1)+sum(lgamma(m+1)))
	if(!log){
		ret <- exp(ret)
	}
	return(ret)
}


#' @export

my.prop.trend.test <- function(m,score=c(0,1,2)){
	x <- m[1:3]
	n <- m[1:3]+m[4:6]
	#tmp <- prop.trend.test(case,all,score)
	freq <- x/n
	p <- sum(x)/sum(n)
	w <- n/p/(1 - p)
	lm.out <- lm(freq ~ score, weights = w)
	return(lm.out[[1]][2])
}

#' @export
my.Cochran.Armitage.Trend.test <- function(m){
	x <- m[1:3]
	n <- m[1:3]+m[4:6]
	X <- sum(x)
	N <- sum(n)
	(N-1)*(N*(x[1]+2*x[2])-X*(n[1]+2*n[2]))^2/(X*(N-X)*(N*(n[1]+4*n[2])-(n[1]+2*n[2])^2))
}

#' @export

my.trend.exact.2x3 <- function(m){
	out1 <- my.all.2x3(m)
	# marginal counts
	sum.col <- c(m[1]+m[4],m[2]+m[5],m[3]+m[6])
	sum.row <- c(sum(m[1:3]),sum(m[4:6]))
	sum.all <- sum(m)

	tmp1 <-sum(lgamma(sum.col+1))+sum(lgamma(sum.row+1))-lgamma(sum.all+1)
	tmp2 <- apply(lgamma(out1+1),1,sum)
	tmp3 <- tmp1 - tmp2
	
	out3 <- apply(out1,1,my.prop.trend.test)
	original.chisq <- my.prop.trend.test(m)
	selected <- which(abs(out3)>=abs(original.chisq))
	sum(exp(tmp3[selected]))
}
