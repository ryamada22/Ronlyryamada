#' 3D triangle mesh smoothing
#'
#' These functions make triangle mesh smooth with Catmull-Clark subdivision
#' @param x.vet is a list. 
#' @keywords shape triangulation mesh smoothing
#' @export
#' @examples
#' library(rgl)
#' library(misc3d)
#' v <- matrix(c(1,2,3,2,3,1,3,3,3),byrow=TRUE,ncol=3)
#' u <- v[c(3,2,1,3,2,2,2,3),]
#' my.tri.v.id(v,u)
#' k <- 4
#' out <- my.3d.obj(3,k)
#' obj <- out$xyz[which(out$xyz.val==1),]
#' plot3d(obj)
#' xyz.arr <- array(out$xyz.val,c(2^k,2^k,2^k))
#' tris <- my.peri.tri(xyz.arr)
#' drawScene.rgl(tris)
#' sss <- my.tri.vet(ttt[[2]])
#' x.vet <- list(ttt[[1]],sss)
#' new.vet <- my.catmull.clark.tri(x.vet)
#' plot3d(new.vet[[1]][new.vet[[2]]$v.of.t,])

my.catmull.clark.tri <- function(x.vet){
	x <- x.vet[[1]]
	v.of.t <- x.vet[[2]]$v.of.t
	t.of.e <- x.vet[[2]]$t.of.e
	v.of.e <- x.vet[[2]]$v.of.e
	face.pt <- (x[v.of.t[,1],] + x[v.of.t[,2],] + x[v.of.t[,3],])/3
	edge.pt <- (face.pt[t.of.e[,1],] + face.pt[t.of.e[,2],] + x[v.of.e[,1],] + x[v.of.e[,2],])/4
	new.v.pt <- matrix(0,length(x[,1]),3)
	for(i in 1:length(new.v.pt[,1])){
		ori.v <- x[i,]
		fs <- x.vet[[2]]$t.of.v[[i]]
		n <- length(fs)
		mean.f.pt <- apply(face.pt[fs,],2,mean)
		ed <- x.vet[[2]]$e.of.v[[i]]
		mean.e.pt <- apply(edge.pt[ed,],2,mean)
		new.v.pt[i,] <- (mean.f.pt + 2*mean.e.pt + (n-3)*ori.v)/n
	}
	new.v.x <- rbind(new.v.pt,face.pt,edge.pt)
	n.v <- length(new.v.pt[,1])
	n.f <- length(face.pt[,1])
	
	new.tris <- matrix(0,n.f*6,3)
	for(i in 1:length(v.of.t[,1])){
		e.id <- x.vet[[2]]$e.of.t[i,]
		v.of.e.of.t <- v.of.e[e.id,]
		
		new.f.id <- i + n.v
		new.e.id <- e.id + n.v + n.f
		
		tmp <- matrix(0,6,3)
		tmp[1:2,1] <- v.of.e.of.t[1,]
		tmp[3:4,1] <- v.of.e.of.t[2,]
		tmp[5:6,1] <- v.of.e.of.t[3,]
		tmp[1:2,2] <- new.e.id[1]
		tmp[3:4,2] <- new.e.id[2]
		tmp[5:6,2] <- new.e.id[3]
		tmp[,3] <- new.f.id
		
		new.tris[(1:6) + (i-1)*6,] <- tmp
	}
	new.tris <- t(apply(new.tris,1,sort))
	new.vet <- my.tri.vet(new.tris)
	return(list(new.v.x,new.vet))
}

#' @export
my.tri.v.id <- function(v,u){
	d <- matrix(0,length(v[,1]),length(u[,1]))
	for(i in 1:length(v[1,])){
		d <- d+abs(outer(v[,i],u[,i],"-"))
	}
	ret <- which(d==0,arr.ind=TRUE)
	ret[,c(2,1)]
}

#' @export
my.tri.vid <- function(tris){
	v <- unique(rbind(tris$v1,tris$v2,tris$v3))
	v1 <- my.tri.v.id(v,tris$v1)
	v2 <- my.tri.v.id(v,tris$v2)
	v3 <- my.tri.v.id(v,tris$v3)
	tri.vid <- cbind(v1[,2],v2[,2],v3[,2])
	tri.vid <- t(apply(tri.vid,1,sort))
	return(list(v=v,tri.vid=tri.vid,tri.coord=list(tris$v1,tris$v2,tris$v3)))
}

#' @export
my.tri.edge <- function(tri){
	ed <- rbind(tri[,1:2],tri[,2:3],tri[,c(1,3)])
	ed <- t(apply(ed,1,sort))
	ed <- unique(ed)
	ed
}

#' @export
my.tri.vet <- function(tri){
	v.of.e <- my.tri.edge(tri)
	n.v <- length(unique(c(tri)))
	#v <- 1:length((x[,1]))
	v.of.t <- t(apply(tri,1,sort))
	e.of.v <- t.of.v <- list()
	e.of.t <- matrix(0,length(tri[,1]),3)
	t.of.e <- matrix(0,length(v.of.e[,1]),2)
	for(i in 1:n.v){
		e.of.v[[i]] <- sort(c(which(v.of.e[,1]==i),which(v.of.e[,2]==i)))
		t.of.v[[i]] <- sort(c(which(v.of.t[,1]==i),which(v.of.t[,2]==i),which(v.of.t[,3]==i)))
	}
	tmp.et.mat <- matrix(0,length(v.of.e[,1]),length(tri[,1]))
	for(i in 1:length(v.of.e[,1])){
		tmp1 <- t(tri[,1:2])-v.of.e[i,]
		tmp2 <- t(tri[,c(1,3)])-v.of.e[i,]
		tmp3 <- t(tri[,2:3])-v.of.e[i,]
		tmp123 <- apply(abs(tmp1),2,sum)*apply(abs(tmp2),2,sum)*apply(abs(tmp3),2,sum)
		tmp4 <- which(tmp123==0)
		tmp.et.mat[i,tmp4] <- 1		
	}
	for(i in 1:length(tmp.et.mat[,1])){
		t.of.e[i,] <- which(tmp.et.mat[i,]==1)
	}
	for(i in 1:length(tmp.et.mat[1,])){
		e.of.t[i,] <- which(tmp.et.mat[,i]==1)
	}
	return(list(e.of.v=e.of.v,t.of.v=t.of.v,v.of.e=v.of.e,t.of.e=t.of.e,v.of.t=v.of.t,e.of.t=e.of.t))
}

