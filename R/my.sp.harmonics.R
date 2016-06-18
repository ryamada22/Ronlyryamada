#' Spherical harmonics
#'
#' Spherical harmonics

#' @param x, a matrix; num. of column = dimension
#' @param R is a real value of radius
#' @export
#' @examples
#' require(rgl)
#' require(RFOC)
#' require(Ronlyryamada)
#' k <- 20
#' Pms <- my.Legendre(k)
#' x <- seq(from=-1,to=1,length=100)
#' ys <- lapply(Pms,my.poly.calc,x)
#' plot(x,ys[[1]],type="l",ylim=c(-1,1))
#' for(i in 2:length(ys)){
#' 	points(x,ys[[i]],type="l",col=i)
#' }
#' n <- 6
#' A <- matrix(0,n,n)
#' ret <- matrix(list(),n,n)
#' for(i in 1:n){
#' 	for(j in 1:i){
#' 		A. <- A
#' 		A.[i,j] <- 1
#' 		ret[i,j] <- list(my.spherical.harm(A.))
#' 	}
#' }
#' tmp <- ret[4,1][[1]]
#' out <- rbind(tmp$X,rep(min(tmp$X),3),rep(max(tmp$X),3))
#' plot3d(out,type="l",col=rainbow(128))
#' col1 <- col2 <- rep(0,length(tmp$r))
#' col1[which(tmp$r>=0)] <- tmp$r[which(tmp$r>=0)]
#' col2[which(tmp$r<0)] <- -tmp$r[which(tmp$r<0)]
#' col1 <- col1/max(col1)
#' col2 <- col2/max(col2)
#' plot3d(tmp$X.,col=rgb(1-col1,1-col2,0),size=10)
#' n <- 6
#' A <- matrix(0,n,n)
#' ret <- matrix(list(),n,n)
#' for(i in 1:n){
#' 	for(j in 1:i){
#' 		A. <- A
#' 		A.[i,j] <- 1
#' 		ret[i,j] <- list(my.spherical.harm(A.))
#' 	}
#' }
#' # 半径そのまま
#' tmp <- ret[4,1][[1]]
#' out <- rbind(tmp$X,rep(min(tmp$X),3),rep(max(tmp$X),3))
#' plot3d(out,type="l",col=rainbow(128))
#' # 半径の２乗で
#' out.. <- rbind(tmp$X*tmp$r,rep(min(tmp$X*tmp$r),3),rep(max(tmp$X*tmp$r),3))
#' plot3d(out..,type="l",col=rainbow(128))
#' # 半径を色で
#' col1 <- col2 <- rep(0,length(tmp$r))
#' col1[which(tmp$r>=0)] <- tmp$r[which(tmp$r>=0)]
#' col2[which(tmp$r<0)] <- -tmp$r[which(tmp$r<0)]
#' col1 <- col1/max(col1)
#' col2 <- col2/max(col2)
#' plot3d(tmp$X.,col=rgb(1-col1,1-col2,0),size=10)
#' n <- 6
#' A <- matrix(0,n,n)
#' ret <- matrix(list(),n,n)
#' for(i in 1:n){
#' 	for(j in 1:i){
#' 		A. <- A
#' 		A.[1,1] <- 5 # 正円の分
#' 		A.[i,j] <- 1
#' 		ret[i,j] <- list(my.spherical.harm(A.))
#' 	}
#' }
#' tmp <- ret[4,1][[1]]
#' out <- rbind(tmp$X,rep(min(tmp$X),3),rep(max(tmp$X),3))
#' plot3d(out,type="l",col=rainbow(128))
#' n <- 6
#' A. <- matrix(runif(n^2),n,n)
#' A.[1,1] <- 5
#' xxx <- my.spherical.harm.mesh(A.,n=32)
#' plot3d(xxx$v)
#' segments3d(xxx$v[c(t(xxx$edge)),])


my.Legendre <- function(k){
	Pm <- list()
	Pm[[1]] <- c(1)
	Pm[[2]] <- c(0,1)
	if(k <= 1){
		return(Pm)
	}
	for(i in 3:(k+1)){
		z <- i-1
		x1 <- c(0,Pm[[i-1]])
		x2 <- Pm[[i-2]]
		L.x1 <- length(x1)
		L.x2 <- length(x2)
		L.max <- max(L.x1,L.x2)
		x1 <- c(x1,rep(0,L.max-L.x1))
		x2 <- c(x2,rep(0,L.max-L.x2))
		Pm[[i]] <- 1/(z)*((2*z-1)*x1 - (z-1)*x2)
	}
	Pm
}

#' @export

my.poly.calc <- function(Pm,x){
	n <- length(Pm)
	ret <- rep(0,length(x))
	for(i in 1:n){
		ret <- ret + Pm[i] * x^(i-1)
	}
	ret
}
#' @export

my.differential <- function(x,k){
	n <- length(x)
	if(k==0){
		ret <- x
	}else if(n < (k+1)){
		ret <- c(0)
	}else{
		tmp <- x[(k+1):n]
		a <- 1:(n-k)
		b <- (k):(n-1)
		ret <- factorial(b)/factorial(a) * tmp
	}
	ret
}
#' @export

my.spherical.harm <- function(A,B=matrix(0,ncol=ncol(A),nrow=nrow(A)),Pm=my.Legendre(length(A[,1])),theta=seq(from=0,to=1,length=100)*2*pi,phi=seq(from=0,to=1,length=100)*pi,normalization=TRUE){
# 角度の座標
	tp <- expand.grid(theta,phi)
	z <- rep(0,length(tp[,1]))
	M <- ncol(A)
# 重みづけ行列でのループ
	for(i in 1:M){
		for(j in 1:i){
			tmp1 <- (A[i,j] * cos((j-1)*tp[,1] + B[i,j]) ) * (-1)^(j-1)*(1-cos(tp[,2])^2)^((j-1)/2)
# ルジャンドル多項式のj-1階微分の係数ベクトル
			tmp2 <- my.differential(Pm[[i]],j-1)
			tmp3 <- cos(tp[,2])
			tmp4 <- rep(0,length(z))
			for(k in 1:length(tmp2)){
				tmp4 <- tmp4 + tmp3^(k-1) * tmp2[k]
			}
			if(normalization){
				tmp5 <- sqrt((2*(i-1)+1)/(4*pi)*factorial(i-j)/factorial(i+j-2))
			}else{
				tmp5 <- 1
			}
			
			z <- z + tmp1 * tmp4 * tmp5
		}
	}
	X <- z * cos(tp[,1]) * sin(tp[,2])
	Y <- z * sin(tp[,1]) * sin(tp[,2])
	Z <- z * cos(tp[,2])
	X. <- cos(tp[,1]) * sin(tp[,2])
	Y. <- sin(tp[,1]) * sin(tp[,2])
	Z. <- cos(tp[,2])
	return(list(r=z,tp=tp,X=cbind(X,Y,Z),X.=cbind(X.,Y.,Z.)))
}

#' @export

my.spherical.harm.mesh <- function(A,B=matrix(0,ncol=ncol(A),nrow=nrow(A)),Pm=my.Legendre(length(A[,1])),n=32,normalization=TRUE){
	tmp <- my_sphere_tri_mesh(n)
	xyz <- tmp$xyz
	tri <- tmp$triangles
	edge <- tmp$edge
	tosp <- TOSPHERE(xyz[,1],xyz[,2],xyz[,3])
	tp <- cbind(tosp[[1]]/360*2*pi,tosp[[2]]/360*2*pi)
	z <- rep(0,length(tp[,1]))
	M <- ncol(A)
# 重みづけ行列でのループ
	for(i in 1:M){
		for(j in 1:i){
			tmp1 <- (A[i,j] * cos((j-1)*tp[,1] + B[i,j]) ) * (-1)^(j-1)*(1-cos(tp[,2])^2)^((j-1)/2)
# ルジャンドル多項式のj-1階微分の係数ベクトル
			tmp2 <- my.differential(Pm[[i]],j-1)
			tmp3 <- cos(tp[,2])
			tmp4 <- rep(0,length(z))
			for(k in 1:length(tmp2)){
				tmp4 <- tmp4 + tmp3^(k-1) * tmp2[k]
			}
			if(normalization){
				tmp5 <- sqrt((2*(i-1)+1)/(4*pi)*factorial(i-j)/factorial(i+j-2))
			}else{
				tmp5 <- 1
			}
			
			z <- z + tmp1 * tmp4 * tmp5
		}
	}
	X <- z * cos(tp[,1]) * sin(tp[,2])
	Y <- z * sin(tp[,1]) * sin(tp[,2])
	Z <- z * cos(tp[,2])
	return(list(v=cbind(X,Y,Z),f=tri,edge=edge))
}



