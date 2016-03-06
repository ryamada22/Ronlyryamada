#' 3D+time movement
#'
#' 3D + time movement
#' @keywords movement
#' @export
#' @export
#' @examples
#' library(rgl)
#' Nobs <- 5;Nt <- 100; d <- 3
#' X = my.random.walk(Nobs,Nt,d)
#' X.k <- my.k.diff(X)
#' X.d <- my.array.dist(X)

my.random.walk <- function(Nobs,Nt,d){
	Tr <- array(0,c(Nobs,Nt,d))
	for(i in 1:Nobs){
		Tr[i,,] <- matrix(rnorm(Nt*d))
		Tr[i,,] <- apply(Tr[i,,],2,cumsum)
		tmp <- apply(Tr[i,,],2,mean)
		Tr[i,,] <- t(t(Tr[i,,]) + tmp)
	}
	Tr
}

#' @export
my.gravity.walk <- function(Nobs,Nt,d,w=rep(1,Nobs),dt=0.01,dt2=0.01,sd=1){
	x <- matrix(rnorm(Nobs*d),ncol=d)
	v <- matrix(rnorm(Nobs*d),ncol=d)
	u <- my.gravity(x,w)*dt2
	ret <- array(0,c(Nobs,Nt+1,d))
	ret[,1,] <- x
	for(i in 1:Nt){
		v <- v + u
		x <- x + v*dt
		u <- my.gravity(x,w)*dt2
		ret[,i+1,] <- x + matrix(rnorm(Nobs*d,sd=sd),ncol=d)
	}
	return(list(x=ret,v=v,u=u))
}

#' @export
my.plot.traj <- function(X){
	plot3d(X[1,,],type="l")
	tmp <- my.series.seg(length(X[1,,1]))
	for(i in 2:dim(X)[1]){
		segments3d(X[i,tmp,],col=i)
	}
}

#' @export
my.gravity <- function(x,w){
	ret <- x
	for(i in 1:length(x[,1])){
		tmp <- t(t(x)-x[i,])
		d <- apply(tmp^2,1,sum)
		d[which(d==0)] <- 1
		tmpw <- w[i] * w/d
		ret[i,] <- apply(tmp * tmpw,2,sum)
	}
	ret
}

#' @export
my.array.dist <- function(X,method="euclidean",a.ob=1,a.t=2){
	dm <- dim(X)
	Nobs <- dm[1]
	Nt <- dm[2]
	ret <- array(0,c(Nt,Nobs,Nobs))
	for(i in 1:Nt){
		tmp <- X[,i,]
		ret[i,,] <- as.matrix(dist(tmp,method=method))
	}
	ret
}	
#' @export
my.k.diff <- function(X,k=2,shrink=TRUE){
	d <- dim(X)
	Nobs <- d[1]
	ret <- array(0,c(k+1,d))
	ret[1,,,] <- X

	for(i in 1:k){
		st <- 1 + ((i+1) %% 2)
		end <- d[2] - (i %% 2)
		for(j in 1:Nobs){
			tmp <- apply(ret[i,j,,],2,diff)
			ret[i+1,j,st:end,] <- tmp
		}
	}
	if(shrink){
		tmp1 <- k %% 2
		tmp2 <- (k-tmp1)/2
		tmp3 <- (k+tmp1)/2
		ret <- ret[,,(1+tmp2):(d[2]-tmp3),]
	}
	ret
}

#' @export
my.series.seg <- function(n){
	tmp <- cbind(1:(n-1),2:n)
	return(c(t(tmp)))
}

