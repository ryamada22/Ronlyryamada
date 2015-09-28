library(geometry)
library(onion)

my.circum.center <- function(v1,v2,v3){
	e1 <- v2-v1
	e2 <- v3-v2
	e3 <- v1-v3
	S <- sqrt(sum(extprod3d(e1,e2)^2))/2
	
	A <- sum(e1^2)
	B <- sum(e2^2)
	C <- sum(e3^2)
	
	(A*(B+C-A)*v3 + B*(C+A-B)*v1 + C*(A+B-C)*v2)/(16*S^2)
}

v1 <- runif(3)
v2 <- runif(3)
v3 <- runif(3)
circm <- my.circum.center(v1,v2,v3)

sum((circm-v1)^2)
sum((circm-v2)^2)
sum((circm-v3)^2)


my.check.cot <- function(v1,v2,v3,v4){
	alpha <- acos(sum((v1-v2)*(v3-v2))/sqrt(sum((v1-v2)^2))/sqrt(sum((v3-v2)^2)))
	beta <- acos(sum((v1-v4)*(v3-v4))/sqrt(sum((v1-v4)^2))/sqrt(sum((v3-v4)^2)))
	ret1 <- (1/tan(alpha) + 1/tan(beta))/2
	print(1/tan(alpha))
	print(1/tan(beta))
	circm1 <- my.circum.center(v1,v2,v3)
	circm2 <- my.circum.center(v1,v3,v4)
	
	L1 <- sqrt(sum((circm1-(v3+v1)/2)^2))*sign(1/tan(alpha))
	L2 <- sqrt(sum((circm2-(v3+v1)/2)^2))*sign(1/tan(beta))
	ret2 <- (L1 + L2)/sqrt(sum((v3-v1)^2))
	#ret2 <- sqrt(sum((circm1-circm2)^2))
	return(list(ret1,ret2))
}

n.iter <- 1000
diffs <- rep(NA,n.iter)
for(i in 1:n.iter){
v1 <- runif(3)
v2 <- runif(3)
v3 <- runif(3)
v4 <- runif(3)

out <- my.check.cot(v1,v2,v3,v4)
diffs[i] <- out[[1]]-out[[2]]
}

my.circ.center.bary <- function(pts){
	d <- length(pts[1,])
	n.pt <- length(pts[,1])
	A <- 2* pts %*% t(pts)
	A <- cbind(A,rep(1,n.pt))
	A <- rbind(A,c(rep(1,n.pt),0))
	b <- c(apply(pts^2,1,sum),1)
	solve(A,b)[-(n.pt+1)]
}

pts <- matrix(c(0,0,4,0),byrow=TRUE,ncol=2)
my.circ.center.bary(pts)

my.circ.center2 <- function(pts){
	tmp <- my.circ.center.bary(pts)
	apply(pts * tmp,2,sum) 
}

my.circ.center2(pts)
pts <- rbind(v1,v2,v3)
my.circum.center(v1,v2,v3)
kkk <- my.circ.center2(pts)
