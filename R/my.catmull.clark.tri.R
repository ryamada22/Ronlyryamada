#' 3D triangle mesh smoothing
#'
#' These functions make triangle mesh smooth with Catmull Clark subdivision
#' @param x.vet is a list. 
#' @keywords shape triangulation mesh smoothing
#' @export
#' @examples
#' library(rgl)
#' library(misc3d)
#' library(Matrix)
#' library(onion)
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
#' ttt <- my.tri.vid(tris)
#' sss <- my.tri.vet(ttt[[2]])
#' x.vet <- list(ttt[[1]],sss)
#' new.vet <- my.catmull.clark.tri(x.vet)
#' plot3d(new.vet[[1]][new.vet[[2]]$v.of.t,])
#' drawScene.rgl(tris)
#' tris2 <- tris
#' tris2$v1 <- new.vet[[1]][new.vet[[2]]$v.of.t[,1],]
#' tris2$v2 <- new.vet[[1]][new.vet[[2]]$v.of.t[,2],]
#' tris2$v3 <- new.vet[[1]][new.vet[[2]]$v.of.t[,3],]
#' drawScene.rgl(tris2)
#' vot <- new.vet[[2]]$v.of.t
#' voe <- new.vet[[2]]$v.of.e
#' toe <- new.vet[[2]]$t.of.e
#' tri.sorted <- my.sort.tri.dir(toe,voe,vot)
#' vertices <- new.vet[[1]][,1]*Hi+new.vet[[1]][,2]*Hj+new.vet[[1]][,3]*Hk
#' faces.v <- t(tri.sorted)
#' my.mesh.tri.plot(vertices,faces.v)
#' tmp.out <- my.mesh.back.euler(vertices,faces.v,step=0.005,eps=10^(-10),max.iter=200)
#' for(i in 1:20){
#' 	tmp.out <- my.mesh.back.euler(tmp.out,faces.v,step=0.005)
#' 	my.mesh.tri.plot(tmp.out,faces.v)
#' }
#' my.mesh.tri.plot(tmp.out,faces.v)


my.catmull.clark.tri <- function(x.vet){
	x <- x.vet[[1]]
	v.of.t <- x.vet[[2]]$v.of.t
	t.of.e <- x.vet[[2]]$t.of.e
	v.of.e <- x.vet[[2]]$v.of.e
	face.pt <- (x[v.of.t[,1],] + x[v.of.t[,2],] + x[v.of.t[,3],])/3
	face.pt <- matrix(face.pt,ncol=3)
	edge.pt <- (face.pt[t.of.e[,1],] + face.pt[t.of.e[,2],] + x[v.of.e[,1],] + x[v.of.e[,2],])/4
	new.v.pt <- matrix(0,length(x[,1]),3)
	for(i in 1:length(new.v.pt[,1])){
		ori.v <- x[i,]
		fs <- x.vet[[2]]$t.of.v[[i]]
		n <- length(fs)
		if(n==1){
			mean.f.pt <- face.pt[fs,]
		}else{
			mean.f.pt <- apply(face.pt[fs,],2,mean)
		}
		
		ed <- x.vet[[2]]$e.of.v[[i]]
		if(length(ed)==1){
			mean.e.pt <- edge.pt[ed,]
		}else{
			mean.e.pt <- apply(edge.pt[ed,],2,mean)
		}
		
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
	if(!is.matrix(u)){
		u <- matrix(u,nrow=1)
	}
	d <- matrix(0,length(v[,1]),length(u[,1]))
	for(i in 1:length(v[1,])){
		d <- d+abs(outer(v[,i],u[,i],"-"))
	}
	ret <- which(d==0,arr.ind=TRUE)
	ret1 <- ret[,c(2,1)]
	if(!is.matrix(ret1)){
		ret1 <- matrix(ret1,ncol=2)
	}
	return(ret1)
}

#' @export
my.tri.vid <- function(tris){
	v <- unique(rbind(tris$v1,tris$v2,tris$v3))
	v1 <- my.tri.v.id(v,matrix(tris$v1,ncol=3))
	v2 <- my.tri.v.id(v,matrix(tris$v2,ncol=3))
	v3 <- my.tri.v.id(v,matrix(tris$v3,ncol=3))
	tri.vid <- cbind(v1[,2],v2[,2],v3[,2])
	tri.vid <- matrix(tri.vid,ncol=3)
	tri.vid.s <- t(apply(tri.vid,1,sort))
	return(list(v=v,tri.vid=tri.vid.s,tri.vid.unsorted=tri.vid,tri.coord=list(tris$v1,tris$v2,tris$v3)))
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
	if(!is.matrix(tri)){
		tri <- matrix(tri,ncol=3)
	}
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
		tmp1 <- t(matrix(tri[,1:2],ncol=2))-v.of.e[i,]
		tmp2 <- t(matrix(tri[,c(1,3)],ncol=2))-v.of.e[i,]
		tmp3 <- t(matrix(tri[,2:3],ncol=2))-v.of.e[i,]
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

#' @export
my.sort.tri.dir <- function(t.of.e,v.of.e,v.of.t){
	n.t <- max(t.of.e)
	n.e <- length(t.of.e[,1])
	n.v <- max(v.of.e)
	tmp <- v.of.e[,1] * n.v + v.of.e[,2]
	ord <- order(tmp)
	mat <- matrix(0,n.t,n.e)
	for(i in 1:n.e){
		mat[t.of.e[i,],i] <- 1
	}
	mat <- mat[,ord]
	ret <- matrix(0,n.t,n.e)
	ret[1,which(mat[1,]==1)] <- c(1,-1,1)

	mat[1,] <- 0
	ori <- rep(0,n.t)
	ori[1] <- 1
	
	loop <- TRUE
	
	while(loop){
		tmp <- apply(ret,2,sum)
		candidate.t <- which(apply(mat,1,sum)==3)
		if(length(candidate.t)==0){
			loop <- FALSE
		}else{
			tmp.ret <- ret
			tmp.ret[candidate.t,] <- -matrix(rep(tmp,length(candidate.t)),byrow=TRUE,ncol=n.e) * mat[candidate.t,]
			diff.ret <- tmp.ret - ret
			abs.diff.ret <- abs(diff.ret)
			new.t <- which(apply(abs.diff.ret,1,sum)>0)
			for(i in 1:length(new.t)){
				eds <- which(mat[new.t[i],]==1)
				tmp.three <- tmp.ret[new.t[i],eds]
				if(prod(tmp.three-c(1,-1,1))==0){
					tmp.ret[new.t[i],eds] <- c(1,-1,1)
					ori[new.t[i]] <- 1
				}else{
					tmp.ret[new.t[i],eds] <- c(-1,1,-1)
					ori[new.t[i]] <- -1
				}
				mat[new.t[i],] <- 0
			}
		}
		ret <- tmp.ret
	}
	new.v.of.t <- v.of.t
	new.v.of.t[which(ori==-1),] <- v.of.t[which(ori==-1),c(1,3,2)]
	
	new.v.of.t
}
#' @export
my.sort.tri.dir.slow <- function(eot,toe,vot,voe){
	n.t <- length(eot[,1])
	n.e <- length(toe[,1])
	ret <- matrix(0,n.t,3)
	done.t <- rep(0,n.t)
	done.e <- rep(0,n.e)
	ret[1,] <- vot[1,]
	seed.t <- 1
	done.t[seed.t] <- 1
	done.e[eot[seed.t,]] <- 1
	loop <- TRUE
	while(loop){
		target.t <- c()
		given.t <- c()
		given.e <- c()
		for(i in 1:length(seed.t)){
			target.e <- eot[seed.t[i],]
			for(j in 1:length(target.e)){
				tmp <- toe[target.e[i],]
				tmp2 <- which(tmp!=seed.t[i])
				target.t <- c(target.t,tmp[tmp2])
				given.t <- c(given.t,seed.t[i])
				given.e <- c(given.e,target.e[j])
			}
		}
		seed.t <- c()
		for(i in 1:length(target.t)){
			if(done.t[i]==0){
				seed.t <- c(seed.t,target.t[i])
				tmp <- ret[given.t[i],]
				ed <- rbind(vot[target.t[i],1:2],vot[target.t[i],2:3],vot[target.t[i],c(1,3)])
				
				ED <- rbind(vot[given.t[i],1:2],vot[given.t[i],2:3],vot[given.t[i],c(1,3)])
				diff1 <- outer(ed[,1],ED[,1],"-")
				diff2 <- outer(ed[,2],ED[,2],"-")
				Diff <- abs(diff1) + abs(diff2)
				id <- which(Diff==0,arr.ind=TRUE)
				ida <- id[1]<3
				idb <- id[2]<3
				if(xor(ida,idb)){
					ret[target.t[i],] <- vot[target.t[i],]
				}else{
					tmp2 <- vot[target.t[i],]
					ret[target.t[i],] <- tmp2[c(2,1,3)]
				}
			}
		}
		done.t[seed.t] <- 1
		done.e[eot[seed.t,]] <- 1

		if(sum(done.e)==length(done.e)){
			loop <- FALSE
		}
	}
	ret
}

#' @export
my.mesh.back.euler <- function(vertices,faces.v,step = 0.01,eps = 10^-5,max.iter = 50){
	# 曲率,rho と法線ベクトル,N
	rho.out <- my.curvature.cot(vertices,faces.v)
	rho.N <- rho.out[[1]]
	new.vertices <- vertices - 2 * rho.N * step
	for(i in 1:max.iter){
		# new.verticesでのrho=new.rho,N=new.N
		rho.out <- my.curvature.cot(new.vertices,faces.v)
		new.rho.N <- rho.out[[1]]
		tmp.vertices <- vertices - 2 * new.rho.N * step
		diff.vertices <- tmp.vertices - new.vertices
		if(max(Norm(diff.vertices) < eps)){
			break
		}
		new.vertices = tmp.vertices		
	}
	new.vertices
}
