#' Utilities for DEC conformal
#'
#' These functions are for DEC conformal

#' @export
rho.fromVtoTri <- function(rho.v,faces.v){
  tmp.rho <- matrix(rho.v[faces.v],nrow=3)
  apply(tmp.rho,2,mean)
}

#' @export
my.vector.access <- function(v,a,func=sum,zero=0){
  if(is.vector(v)){
		v <- matrix(v,ncol=1)
	}
	ord <- order(a)
	rle.out <- rle(a[ord])
	num.row <- length(rle.out[[1]])
	num.col <- max(rle.out[[1]])
	tmp1 <- rep(1:num.row,rle.out[[1]])
	tmp2 <- c()
	for(i in 1:num.row){
		tmp2 <- c(tmp2,1:rle.out[[1]][i])
	}
	addr <- tmp1 + num.row*(tmp2-1)
	ret.v <- matrix(0,num.row,ncol(v))
	for(i in 1:ncol(v)){
		if(zero==0){
			tmp.V <- sparseVector(v[ord,i],i=addr,length=num.row*num.col)
			M <- Matrix(tmp.V,num.row,num.col)
		}else{
			M <- matrix(zero,num.row,num.col)
			M[addr] <- v[ord,i]

		}
		ret.v[,i] <- apply(M,1,func)

	}
	return(list(A = rle.out[[2]],V = ret.v))
}
#' Making E
#' @export
my.make.E.v <- function(vertices,faces.v,rho){
	edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
	edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	tmp <- edge1 * edge2
  A <- Mod(Im(tmp))/2
	coef.a <- -1/(4*A)
	coef.b <- rho/6
	coef.c <- A*rho^2/9
	
	E.re <- E.i <- E.j <- E.k <- sparseVector(c(0),i=c(1),length=length(vertices)^2)
  
	e.q <- list()
	e.q[[1]] <- vertices[faces.v[2,]]-vertices[faces.v[3,]]
	e.q[[2]] <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	e.q[[3]] <- vertices[faces.v[1,]]-vertices[faces.v[2,]]
	for(i in 1:3){
		for(j in 1:3){
			tmp <- coef.a * e.q[[i]] * e.q[[j]]+ coef.b * (e.q[[j]] -e.q[[i]] ) + coef.c
			addr <- faces.v[i,] + (length(vertices)*(faces.v[j,]-1))
			
			tmp.v <- t(as.matrix(tmp))
			tmp.out <- my.vector.access(tmp.v,addr)
			E.re <- E.re + sparseVector(tmp.out[[2]][,1],tmp.out[[1]],length(vertices)^2)
			E.i <- E.i + sparseVector(tmp.out[[2]][,2],tmp.out[[1]],length(vertices)^2)
			E.j <- E.j + sparseVector(tmp.out[[2]][,3],tmp.out[[1]],length(vertices)^2)
			E.k <- E.k + sparseVector(tmp.out[[2]][,4],tmp.out[[1]],length(vertices)^2)
		}
	}

	return(list(E.re=E.re,E.i=E.i,E.j=E.j,E.k=E.k))
}

#' @export
my.qMtorM <- function(Es){
  n <- sqrt(length(Es[[1]]))
	N <- (n*4)^2
	init.id <- c(1:4,(1:4)+n*4,(1:4)+n*4*2,(1:4)+n*4*3)
	spacing.id <- c(outer((0:(n-1)*4),n*4*4*(0:(n-1)),"+"))
	ret <- sparseVector(c(0),i=c(1),N)
	a <- c(1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1)
	b <- c(1,1,1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1)
	for(j in 1:length(a)){
		tmp.v <- sparseVector(b[j] * Es[[a[j]]]@x,i=init.id[j]+spacing.id[Es[[a[j]]]@i],length=N)
		ret <- ret + tmp.v
	}
	Matrix(ret,n*4,n*4)
}
#' @export
my.inv.pow.2 <- function(A,n.iter=3,b=rep(1,ncol(A)),log=FALSE){
  x <- b
	if(log){
		x.log <- matrix(0,n.iter+1,ncol(A))
		x.log[1,] <- x
	}
	#x <- x/sqrt(sum(x^2))
	#A. <- solve(A)
	for(i in 1:n.iter){
		x2 <- solve(A,x)
		x <- x2/sqrt(sum(x2^2))
		if(log){
			x.log[i+1,] <- x
		}
		
	}
	if(log){
		return(list(x=x,x.log=x.log))
	}else{
		return(list(x=x,x.log=matrix(0,0,ncol(A))))
	}
}
#' @export
my.make.L <- function(vertices,faces.v){
  n.v <- length(vertices)
	L <- sparseVector(c(0),i=c(1),length=n.v^2)
	for(i in 1:3){
		v.ord <- ((1:3)+i+1) %% 3 + 1
		k1 <- faces.v[v.ord[1],]
		k2 <- faces.v[v.ord[2],]
		k3 <- faces.v[v.ord[3],]
		v1 <- vertices[k1]
		v2 <- vertices[k2]
		v3 <- vertices[k3]
		
		u1 <- v2-v1
		u2 <- v3-v1
		u12 <- u1 * u2
		cotAlpha <- (-Re(u12))/Mod(Im(u12))
		addrk2k2 <- k2 + (k2-1)*n.v
		addrk3k3 <- k3 + (k3-1)*n.v
		addrk2k3 <- k2 + (k3-1)*n.v
		addrk3k2 <- k3 + (k2-1)*n.v
		
		addr <- c(addrk2k2,addrk3k3,addrk2k3,addrk3k2)
		
		val <- c(cotAlpha,cotAlpha,-cotAlpha,-cotAlpha)/2
		
		tmp.out <- my.vector.access(val,addr)
		L <- L + sparseVector(tmp.out[[2]][,1],tmp.out[[1]],n.v^2)
	}
	L
}
#' @export
my.make.quatList <- function(L){
  L.q <- list()
  L.q[[1]] <- L
  for(i in 2:4){
    L.q[[i]] <- L*0
  }
  L.q
}
#' @export
my.make.omega <- function(vertices,faces.v,lambda){
  n.v <- length(vertices)
	omega <- rep(0*Hi,n.v)
	for(i in 1:3){
		v.ord <- ((1:3)+i+1) %% 3 + 1
		k1 <- faces.v[v.ord[1],]
		k23 <- rbind(faces.v[v.ord[2],],faces.v[v.ord[3],])
		k23 <- apply(k23,2,sort)
		k2 <- k23[2,]
		k3 <- k23[1,]
		v1 <- vertices[k1]
		v2 <- vertices[k2]
		v3 <- vertices[k3]
		
		edge <- v3-v2
		
		lambda.mat <- matrix(lambda,nrow=4)
		lambda.q <- as.quaternion(lambda.mat)
		lambda.mat.2 <- lambda.mat
		lambda.mat.2[2:4,] <- -lambda.mat.2[2:4,]
		lambda.q. <- as.quaternion(lambda.mat.2)
		lambda1 <- lambda.q[k2]
		lambda1. <- lambda.q.[k2]
		lambda2 <- lambda.q[k3]
		lambda2. <- lambda.q.[k3]
		val <- 1/3 * lambda1. * edge * lambda1 + 1/6 * lambda1. * edge * lambda2 + 1/6 * lambda2. * edge * lambda1 + 1/3 * lambda2. * edge * lambda2

		u1 <- v2-v1
		u2 <- v3-v1
		u12 <- u1 * u2
		cotAlpha <- (-Re(u12))/Mod(Im(u12))

		Val <- cotAlpha * val /2
		Val.m <- t(as.matrix(Val))
		tmp.out2 <- my.vector.access(-Val.m,k2)
		tmp.out3 <- my.vector.access(Val.m,k3)
		omega[tmp.out2[[1]]] <- omega[tmp.out2[[1]]] + as.quaternion(t(tmp.out2[[2]]))
		omega[tmp.out3[[1]]] <- omega[tmp.out3[[1]]] + as.quaternion(t(tmp.out3[[2]]))
	}
	omega.re <- as.matrix(omega)
	omega.re <- omega.re-apply(omega.re,1,mean)
	c(omega.re)
}
#' @export
plot.sp.conformal <- function(xyz,faces.v,sp.tri,rho.f,col1=c(4,5)){
  plot3d(xyz,xlab="x",ylab="y",zlab="z")
  mesh.tri <- tmesh3d(t(xyz),faces.v,homogeneous=FALSE)
  
  col. <- rep(col1,length(sp.tri))[1:length(sp.tri)]
  col <- rep(col.,sp.tri*3)
  rho.f <- rep(rho.f,each=3)
  rho.f <- (rho.f-min(rho.f))/(max(rho.f)-min(rho.f))
  
  col2 <- rgb(1-rho.f,1,col/6)
  #col3 <- gray((rho.f2+1)*0.5)
  shade3d(mesh.tri,col=col2)  
}

#' @export
my.conformal.rho <- function(vertices,faces.v,rho.v,face=FALSE){
  xyz.ori <- t(as.matrix(vertices)[2:4,])
  n.v <- length(vertices) 
  n.f <- length(faces.v[1,])
  if(!face){
    rho.f <- rho.fromVtoTri(rho.v,faces.v)
  }else{
    rho.f <- rho.v
  }
  
  
  edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
	edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	tmp <- edge1 * edge2
  A <- Mod(Im(tmp))/2
  
  s <- sum(A*rho.f)
  #rho.f <- rho.f -s*A/sum(A)
  rho.f <- rho.f -s/A/length(A)
  
  E <- my.make.E.v(vertices,faces.v,rho.f)
  E.re <- my.qMtorM(E)
  
  lambda.v <- my.inv.pow.2(E.re)[[1]]
  
  L <- my.make.L(vertices,faces.v)
  L.q <- my.make.quatList(L)
  L.re <- my.qMtorM(L.q)
  
  omega <- my.make.omega(vertices,faces.v,lambda.v)
  
  new.vertices <- as.quaternion(matrix(solve(L.re,omega),nrow=4))
  xyz.new <- t(as.matrix(new.vertices)[2:4,])
  mean.new <- apply(xyz.new,2,mean)
  xyz.new. <- t(t(xyz.new)-mean.new)
  max.new. <- max(abs(xyz.new.))
  xyz.new.st <- xyz.new./max.new.
  
  new.q <- as.quaternion(t(cbind(rep(0,n.v),xyz.new.st)))
  #ret <- xyz.new.st[,1]*Hi + xyz.new.st[,2]*Hj + xyz.new.st[,3]*Hk
  
  ret <- list(xyz.new=xyz.new.st,xyz.ori=xyz.ori,xyz.new.q=new.q,xyz.ori.q=vertices,faces.v=faces.v,E=E,L=L,lambda.v=lambda.v,omega=omega,rho.f=rho.f)
  ret
}

#' @export
my.sphereConformal <- function(n.psi,rho.fx){
  sp.mesh <- my.sphere.tri.mesh(n.psi)
  vertices.mat <- t(sp.mesh$xyz) 
  vertices <- vertices.mat[1,]*Hi + vertices.mat[2,]*Hj + vertices.mat[3,]*Hk
  edges <- sp.mesh$edge
  faces.v <- t(sp.mesh$triangles)
  rho.v <- rho.fx(sp.mesh$xyz) 
  out <- my.conformal.rho(vertices,faces.v,rho.v)
  ret <- list(xyz.new=out$xyz.new,xyz.ori=sp.mesh$xyz,xyz.new.q=out$xyz.new.q,xyz.ori.q=vertices,faces.v=faces.v,E=out$E,L=out$L,lambda.v=out$lambda.v,omega=out$omega,n.psi=n.psi,rho.fx=rho.fx,rho.v=rho.v,rho.f=out$rho.f,sp.mesh=sp.mesh)
  ret
}
#' @export
my.deform.k.serial <- function(vertices,faces.v,k,n.iter=10){
  v.list <- rho.cot.list <- list()
  v.list[[1]] <- vertices
  rho.cot.list[[1]] <- my.curvature.cot(v.list[[1]],faces.v)
  for(i in 1:n.iter){
    
    tmp.rho.f <- rho.cot.list[[i]][[3]] * Mod(Im(rho.cot.list[[i]][[2]]))
    tmp.out <- my.conformal.rho(v.list[[i]],faces.v,k*tmp.rho.f,face=TRUE)
    v.list[[i+1]] <- tmp.out$xyz.new.q
    rho.cot.list[[i+1]] <- my.curvature.cot(v.list[[i+1]],faces.v)
  }
  return(list(v.list=v.list,rho.cot.list=rho.cot.list,k=k))
}

#' @export
my.deform.k.multi <- function(vertices,faces.v,ks){
  v.list <- rho.cot.list <- list()
  v.list[[1]] <- vertices
  rho.cot.list[[1]] <- my.curvature.cot(v.list[[1]],faces.v)
  for(i in 1:length(ks)){
    
    tmp.rho.f <- rho.cot.list[[1]][[3]] * Mod(Im(rho.cot.list[[1]][[2]]))
    tmp.out <- my.conformal.rho(v.list[[1]],faces.v,ks[i]*tmp.rho.f,face=TRUE)
    v.list[[i+1]] <- tmp.out$xyz.new.q
    rho.cot.list[[i+1]] <- my.curvature.cot(v.list[[i+1]],faces.v)
  }
  return(list(v.list=v.list,rho.cot.list=rho.cot.list,k=k))
}

#' @export
my.mesh.tri.plot <- function(vertices,faces.v,rho.f=NULL){
	if(is.null(rho.f)){
		rho.cot <- my.curvature.cot(vertices,faces.v)
		rho.f <- rho.cot[[3]] * Mod(Im(rho.cot[[2]]))
	}
	xyz <- as.matrix(vertices)[2:4,]
	plot3d(t(xyz))
	mesh.tri <- tmesh3d(xyz,faces.v,homogeneous=FALSE)
	rho.f <- rho.cot[[3]] * Mod(Im(rho.cot[[2]]))
	rho.f <- rep(rho.f,each=3)
	rho.f1 <- rho.f2 <- rep(0,length(rho.f))
	rho.f1[which(rho.f>0)] <- rho.f[which(rho.f>0)]
	rho.f2[which(rho.f<0)] <- -rho.f[which(rho.f<0)]
	#col2 <- rgb(rho.f1/max(rho.f1),rho.f2/max(rho.f2),0.5)
	red <- rep(0,length(rho.f))
	green <- rep(0,length(rho.f))
	if(max(rho.f1)!=0){
		red <- rho.f1/max(rho.f1)
	}
	if(max(rho.f2)!=0){
		green <- rho.f2/max(rho.f2)
	}
	col2 <- rgb(red,green,0.5)
	shade3d(mesh.tri,col=col2)
}

#' @export
my.Laplacian <- function(vertices, faces.v){
	n.v <- length(vertices)
	n.f <- length(faces.v[1,])
	n.e <- n.f * 3/2
	
	edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
	edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	tmp <- edge1 * edge2
	tmp2 <- Mod(Im(tmp))
	A.f <- (tmp2)/2
	val <- rep(A.f,each=3)
	addr <- c(faces.v)
	tmp.out <- my.vector.access(val,addr)
	A.v <- tmp.out[[2]][,1] 
	A.v <- A.v/3
	
}

#' @export
my.curvature.cot <- function(vs,faces.v){
	n.v <- length(vs)
	ret <- rep(0*Hk,n.v)
	edge1 <- vs[faces.v[2,]]-vs[faces.v[1,]]
	edge2 <- vs[faces.v[3,]]-vs[faces.v[1,]]
	tmp <- edge1 * edge2
	tmp2 <- Mod(Im(tmp))
	#tmp2 <- i(tmp)+j(tmp)+k(tmp)
	#inv.f <- which(tmp2<0)
	#faces.v[2:3,inv.f] <- rbind(faces.v[3,inv.f],faces.v[2,inv.f])
	A.f <- (tmp2)/2
	val <- rep(A.f,each=3)
	addr <- c(faces.v)
	tmp.out <- my.vector.access(val,addr)
	A.v <- tmp.out[[2]][,1]

	for(i in 1:3){
		v.ord <- ((1:3)+i+1) %% 3 + 1
		k1 <- faces.v[v.ord[1],]
		k2 <- faces.v[v.ord[2],]
		k3 <- faces.v[v.ord[3],]
		v1 <- vs[k1]
		v2 <- vs[k2]
		v3 <- vs[k3]
		
		u1 <- v2-v1
		u2 <- v3-v1

		u3 <- v3-v2
		#u3.len <- Mod(Im(u3))
		u12 <- u1 * u2
		cotAlpha <- ((-Re(u12))/Mod(Im(u12)))

		val <- c(cotAlpha * u3, cotAlpha * (-u3))
		addr <- c(k2,k3)
		tmp.out <- my.vector.access(t(as.matrix(val)),addr)
		tmp.val <- as.quaternion(t(tmp.out[[2]]))
		ret[tmp.out[[1]]] <- ret[tmp.out[[1]]] + tmp.val
	}
	ret.vec <- -ret/(4*A.v)
	ret.face.re <- rho.fromVtoTri(Re(ret.vec),faces.v)
	ret.face.i <- rho.fromVtoTri(i(ret.vec),faces.v)
	ret.face.j <- rho.fromVtoTri(j(ret.vec),faces.v)
	ret.face.k <- rho.fromVtoTri(k(ret.vec),faces.v)
	ret.face <- ret.face.re + Hi * ret.face.i + Hj * ret.face.j + Hk * ret.face.k


	dir <- sign(Re(tmp * ret.face))
	return(list(norm.v=ret.vec,norm.face=ret.face,dir=dir))
}


