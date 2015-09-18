#' 2D volume and circumference
#'
#' These functions for 2D 0/1 image matrices
#' @param x is a 0/1 matrix.
#' @param r is a positive real value of width to measure circumference.
#' @keywords shape circumference
#' @export
#' @examples
#' x <- matrix(1:16,4,4)
#' v <- c(2,3)
#' 
#' my.parallel(x,v)
#' v <- c(-2,-3)
#' my.parallel(x,v)
#' v <- c(2,-3)
#' my.parallel(x,v)
#' n <- 2^8
#' x <- y <- 1:n
#' xy <- as.matrix(expand.grid(x,y))
#' xy.val <- rep(0,length(xy[,1]))
#' n.circle <- 100
#' ctr <- matrix(sample(50:200,replace=TRUE,n.circle*2),ncol=2)
#' r <- sample(1:40,replace=TRUE,n.circle)
#' 
#' for(i in 1:n.circle){
#' 	tmp.x <- xy[,1] - ctr[i,1]
#' 	tmp.y <- xy[,2] - ctr[i,2]
#' 	s <- which(tmp.x^2+tmp.y^2 < r[i]^2)
#' 	xy.val[s] <- 1
#' }
#' xy.mat <- matrix(xy.val,ncol=n)
#' image(xy.mat)
#' rs <- seq(from=0,to=20,length=21)
#' rs <- rs[-1]
#' tmp.out <- list()
#' for(i in 1:length(rs)){
#' 	tmp.out[[i]] <- my.circ.measure(xy.mat,rs[i])
#' }
#' circs <- rep(0,length(rs))
#' for(i in 1:length(tmp.out)){
#' 	circs[i] <- tmp.out[[i]]$circ
#' }
#' plot(rs,circs)

my.circ.measure <- function(x,r){
	out <- my.sphere.circ(r)
	tmp.x.list <- array(0,c(dim(x),length(out[,1])))
	for(i in 1:length(out[,1])){
		tmp.v <- out[i,]
		tmp.x.list[,,i] <- my.parallel(x,tmp.v)
		
	}

	tt <- apply(tmp.x.list,c(1,2),sum) - length(out[,1])*x
	a <- length(which(tt!=0))
	return(list(band=tt,area=a,width=2*r,circ=a/(2*r)))
}


#' @export
# 格子グラフ


# 原点を中心に半径r以内にある点が作るグラフの外縁の点を列挙
my.sphere.circ <- function(r){
	t <- 0:r
	tmp <- sqrt(r^2-t^2)
	ret <- floor(tmp)
	ret <-cbind(t,ret)
	ret <- rbind(ret,ret[,2:1])
	ret <- unique(ret)
	ret <- rbind(ret,cbind(ret[,1],-ret[,2]),cbind(-ret[,1],ret[,2]),-ret)
	unique(ret)
}




# 平行移動差分
#' @export
my.parallel <- function(x,v){
	dim.x <- dim(x)
	ret <- matrix(0,dim.x[1],dim.x[2])
	if(all(dim.x-abs(v)>0)){
		# コピーするべき範囲
		x.copy <- c(1,dim.x[1]-abs(v[1]))
		y.copy <- c(1,dim.x[2]-abs(v[2]))
		# ペーストするべき範囲
		x.paste <- c(abs(v[1])+1,dim.x[1])
		y.paste <- c(abs(v[2])+1,dim.x[2])
		if(v[1] < 0){
			tmp <- x.copy
			x.copy <- x.paste
			x.paste <- tmp
		}
		if(v[2]<0){
			tmp <- y.copy
			y.copy <- y.paste
			y.paste <- tmp
		}
		ret[x.paste[1]:x.paste[2],y.paste[1]:y.paste[2]] <- x[x.copy[1]:x.copy[2],y.copy[1]:y.copy[2]]
	}else{
		warning("out of range")
	}
	ret
}

