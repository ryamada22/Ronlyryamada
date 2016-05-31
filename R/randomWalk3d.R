#' 3D+time movement
#'
#' 3D + time movement
#' @keywords movement
#' @export
#' @export
#' @examples
#' library(rgl)
#' library(onion)
#' Nobs <- 5;Nt <- 100; d <- 3
#' X = my.random.walk(Nobs,Nt,d)
#' X.k <- my.k.diff(X)
#' X.d <- my.array.dist(X)
#' X.mf <- array(c(t(apply(X.k,2,my.moving.frame.array))),c(Nobs,Nt-2,3,3))



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
my.gravity.walk <- function(Nobs,Nt,d,w=rep(1,Nobs),dt=0.01,dt2=0.01,sd=0.01){
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
	if(shrink){
		ret <- array(0,c(k+1,d))
		ret[1,,,] <- X
		for(i in 1:k){
			st <- 1
			end <- d[2] - 1
			for(j in 1:Nobs){
				tmp <- apply(ret[i,j,,],2,diff)
				ret[i+1,j,st:end,] <- tmp
			}
		}
		ret <- ret[,,1:(d[2]-k),]
		return(ret)
	}
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
	return(ret)
}

#' @export
my.gramSchmidt.mf <- function(X){
	d <- dim(X)
	X.k <- my.k.diff(X,k=d[3])
	Q <- R <- array(0,c(d[1],length(X.k[1,1,,1]),d[3],d[3]))
	for(i in 1:d[1]){
		for(j in 1:length(X.k[1,1,,1])){
			tmp <- gramSchmidt(t(X.k[2:(d[3]+1),i,j,]))
			Q[i,j,,] <- tmp$Q
			R[i,j,,] <- tmp$R
		}
	}
	return(list(Q=Q,R=R,X.k=X.k))
}

#' @export
my.series.seg <- function(n){
	tmp <- cbind(1:(n-1),2:n)
	return(c(t(tmp)))
}

#' @export
my.moving.frame <- function(X,X.,X..){
	H <- my.quaternion.conversion(X)
	H. <- my.quaternion.conversion(X.)
	H.. <- my.quaternion.conversion(X..)
	H... <- H. * H..
	
	H..2 <- H... * H.
	MF1 <- H./Mod(H.)
	MF2 <- H..2/Mod(H..2)
	MF3 <- H.../Mod(H...)
	MF1. <- as.matrix(Im(MF1))[2:4,]
	MF2. <- as.matrix(Im(MF2))[2:4,]
	MF3. <- as.matrix(Im(MF3))[2:4,]
	tmp <- cbind(t(MF1.),t(MF2.),t(MF3.))
	return(array(c(tmp),c(length(MF1),3,3)))
	#return(cbind(MF1,MF2,MF3))
}

#' @export
my.moving.frame.array <- function(X){
	my.moving.frame(X[1,,],X[2,,],X[3,,])
}

#' @export
my.quaternion.conversion <- function(x){
	if(!is.matrix(x)){
		x <- matrix(x,ncol=3)
	}
	Hi * x[,1] + Hj * x[,2] + Hk * x[,3]
}
