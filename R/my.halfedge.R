#' quaternion matrix
#' @export
#' @examples
#' library(Matrix)
#' example(my.catmull.clark.tri)
#' mesh <- my.mesh(faces.v)
#' xyz <- t(as.matrix(vertices)[2:4,])
#' my.mesh.tri.plot(vertices,faces.v)
#' area.norm <- my.area.norm(mesh,xyz)
#' cot <- my.halfedge.cot(mesh,xyz)
#' dec <- my.dec(mesh,xyz)
#' n.iter <- 100
#' current.xyz <- xyz
#' for(i in 1:n.iter){
#' 	new.out <- my.update.heatflow(mesh,current.xyz,step=0.05)
#' 	my.mesh.tri.plot.xyz(new.out$xyz,faces.v)
#' 	current.xyz <- new.out$xyz
#' }



my.update.heatflow <- function(mesh,xyz,step=0.01){
	dec <- my.dec(mesh,xyz)
	L <- t(dec$d0) %*% dec$star1 %*% dec$d0
	#L <- L - (max(L)+min(L))/2
	A <- dec$star0 + step * L
	rhs <- dec$star0 %*% xyz
	xyz <- Matrix::solve(A,rhs)
	xyz <- my.normalize.xyz(xyz)
	dec <- my.dec(mesh,xyz)
	return(list(mesh=mesh,xyz=xyz,dec=dec))
}
#' @export
my.mesh.tri.plot.xyz <- function(xyz,faces.v){
	vertices <- xyz[,1]*Hi+xyz[,2]*Hj + xyz[,3]*Hk
	my.mesh.tri.plot(vertices,faces.v)
}

#' @export
my.normalize.xyz <- function(xyz){
	m <- apply(xyz,2,mean)
	ret <- t(t(xyz)-m)
	M <- sqrt(apply(ret^2,1,sum))
	ret <- ret/max(M)
	ret
}
#' @export

my.qM <- function(m.list){
	d <- dim(m.list[[1]])
	M <- Matrix(0,d[1]*4,d[2]*4)
	tmp <- m.list[[1]]^2 + m.list[[2]]^2 + m.list[[3]]^2 + m.list[[4]]^2
	arr <- which(tmp!=0,arr.ind=TRUE)
	for(i in 1:length(arr[,1])){
		add1 <- (arr[i,1]-1)*4 + 1:4
		add2 <- (arr[i,2]-1)*4 + 1:4
		tmp <- matrix(0,4,4)
		diag(tmp) <- m.list[[1]][arr[i,1],arr[i,2]]
		tmp[2,1] <- tmp[4,3] <- m.list[[2]][arr[i,1],arr[i,2]]
		tmp[1,2] <- tmp[3,4] <- - m.list[[2]][arr[i,1],arr[i,2]]
		tmp[3,1] <- tmp[2,4] <- m.list[[3]][arr[i,1],arr[i,2]]
		tmp[1,3] <- tmp[4,2] <- -m.list[[3]][arr[i,1],arr[i,2]]
		tmp[4,1] <- tmp[3,2] <- m.list[[4]][arr[i,1],arr[i,2]]
		tmp[1,4] <- tmp[2,3] <- -m.list[[4]][arr[i,1],arr[i,2]]

		M[add1,add2] <- tmp
	}
	M
}
#' @export

my.mesh <- function(faces.v){
	n.v <- max(c(faces.v))
	n.f <- length(faces.v[1,])
	n.e <- n.f*3/2
	n.he <- n.e * 2
	HE <- list(id=1:n.he,v=rep(NA,n.he),f=rep(NA,n.he),e=rep(NA,n.he),flip=rep(NA,n.he),nxt=rep(NA,n.he))
	
	F <- list(id=1:n.f,he=matrix(NA,n.f,3),v=matrix(NA,n.f,3))
	V <- list(id=1:n.v,he=rep(NA,n.v),f=list())
	E <- list(id=1:n.e,he=rep(NA,n.e))

	for(i in 1:n.f){
		cnt1 <- (i-1) * 3 + 1
		cnt2 <- (i-1) * 3 + 2
		cnt3 <- (i-1) * 3 + 3
		HE$v[cnt1] <- faces.v[1,i]
		HE$v[cnt2] <- faces.v[2,i]
		HE$v[cnt3] <- faces.v[3,i]
		HE$f[c(cnt1,cnt2,cnt3)] <- i
		HE$nxt[cnt1] <- cnt2
		HE$nxt[cnt2] <- cnt3
		HE$nxt[cnt3] <- cnt1
		F$he[i,] <- c(cnt1,cnt2,cnt3)
		if(is.na(V$he[faces.v[1,i]])){
			V$he[faces.v[1,i]] <- cnt1
		}
		if(is.na(V$he[faces.v[2,i]])){
			V$he[faces.v[2,i]] <- cnt2
		}
		if(is.na(V$he[faces.v[3,i]])){
			V$he[faces.v[3,i]] <- cnt3
		}
	}
	
	v.st <- v.end <- rep(NA,n.he)
	for(i in 1:n.he){
		v.st[i] <- HE$v[i]
		v.end[i] <- HE$v[HE$nxt[i]]
	}
	v.st.end <- t(apply(cbind(v.st,v.end),1,sort))
	v.st.end.value <- v.st.end[,1]-1 + v.st.end[,2]*(n.v-1)
	ord <- order(v.st.end.value)
	ord2 <- matrix(ord,byrow=TRUE,ncol=2)
	HE$flip[ord2[,1]] <- ord2[,2]
	HE$flip[ord2[,2]] <- ord2[,1]
	HE$e[ord2[,1]] <- 1:n.e
	HE$e[ord2[,2]] <- 1:n.e
	E$he <- ord2[,1]

	V$f <- my.enumerate.v.f2(V,HE)
	F$v <- my.enumerate.f.v2(F,HE)

	return(list(HE=HE,V=V,E=E,F=F,n.f=n.f,n.v=n.v,n.e=n.e,n.he=n.he))
}
#' @export

my.area.norm <- function(mesh,xyz){
	ret.face.area <-	rep(NA,mesh$n.f)
	ret.face.norm <- matrix(NA,mesh$n.f,3)
	for(i in 1:mesh$n.f){
		tmp <- my.f.area.norm(xyz[mesh$F$v[i,],])
		ret.face.area[i] <- tmp$area
		ret.face.norm[i,] <- tmp$norm
	}
	ret.vertex.area <- rep(NA,mesh$n.v)
	ret.vertex.norm <- matrix(NA,mesh$n.v,3)
	for(i in 1:mesh$n.v){
		ret.vertex.area[i] <- sum(ret.face.area[mesh$V$f[[i]]])/3
		tmp <- apply(ret.face.norm[mesh$V$f[[i]],],2,sum)
		ret.vertex.norm[i,] <- tmp/sqrt(sum(tmp^2))
	}
	e.he <- mesh$E$he
	e.he.nxt <- mesh$HE$nxt[e.he]
	ret.edge <- xyz[mesh$HE$v[e.he.nxt],] - xyz[mesh$HE$v[e.he],]
	
	return(list(face.area=ret.face.area,face.norm=ret.face.norm,vertex.area=ret.vertex.area,vertex.norm=ret.vertex.norm,edge.vector=ret.edge))
}
#' @export

my.halfedge.cot <- function(mesh,xyz){
	cot <- rep(NA,mesh$n.he)
	for(i in 1:mesh$n.he){
		p0 <- mesh$HE$v[mesh$HE$nxt[mesh$HE$nxt[i]]]
		p1 <- mesh$HE$v[i]
		p2 <- mesh$HE$v[mesh$HE$nxt[i]]
		v1 <- xyz[p1,]-xyz[p0,]
		v2 <- xyz[p2,]-xyz[p0,]
		tmp <- c(v1[2]*v2[3]-v1[3]*v2[2],v1[3]*v2[1]-v1[1]*v2[3],v1[1]*v2[2]-v1[2]*v2[1])
		cot[i] <- sum(v1*v2)/sqrt(sum(tmp^2))
	}
	cot
}
#' @export

my.dec <- function(mesh,xyz,area.norm=my.area.norm(mesh,xyz),cot=my.halfedge.cot(mesh,xyz)){
	star0 <- Diagonal(x=area.norm$vertex.area)
	e.cot <- rep(NA,mesh$n.e)
	d0 <- Matrix(0,mesh$n.e,mesh$n.v)
	for(i in 1:mesh$n.e){
		he <- mesh$E$he[i]
		flip.he <- mesh$HE$flip[he]
		e.cot[i] <- (cot[he] + cot[flip.he])/2
		d0[i,mesh$HE$v[he]] <- -1
		d0[i,mesh$HE$v[flip.he]] <- 1
	}
	star1 <- Diagonal(x=e.cot)
	star2 <- Diagonal(x=area.norm$face.area)
	
	return(list(d0=d0,star0=star0,star1=star1,e.cot=e.cot))
}
#' @export

my.Laplacian <- function(d0,star0,star1){
	solve(star0) %*% t(d0) %*% star1 %*% d0
}
#' @export

my.Laplacian2 <- function(d0,star0,star1){
	t(d0) %*% star1 %*% d0
}
#' @export

my.div <- function(d0,star1){
	t(d0) %*% star1
}
#' @export

my.Dirac <- function(mesh,xyz,area.norm){
	ret <- Matrix(0,mesh$n.f*4,mesh$n.v*4)
	m.list <- list()
	for(i in 1:4){
		m.list[[i]] <- Matrix(0,mesh$n.f,mesh$n.v)
	}
	for(i in 1:mesh$n.f){
		tmp.he <- mesh$F$he[i,]
		tmp.he.nxt <- mesh$HE$nxt[tmp.he]
		tmp.he.nxt2 <- mesh$HE$nxt[tmp.he.nxt]
		tmp.v <- mesh$HE$v[tmp.he.nxt2]
		tmp.he.vector <- xyz[mesh$HE$v[tmp.he.nxt],] - xyz[mesh$HE$v[tmp.he],]
		tmp.area <- area.norm$face.area[i]
		m.list[[2]][i,tmp.v] <- -tmp.he.vector[,1]/tmp.area/2
		m.list[[3]][i,tmp.v] <- -tmp.he.vector[,2]/tmp.area/2
		m.list[[4]][i,tmp.v] <- -tmp.he.vector[,3]/tmp.area/2
	}
	my.qM(m.list)
}
#' @export

my.B <- function(mesh){
	m.list <- list()
	for(i in 1:4){
		m.list[[i]] <- Matrix(0,mesh$n.f,mesh$n.v)
		if(i==1){
			for(j in 1:mesh$n.f){
				m.list[[i]][j,mesh$F$v[j,]] <- 1/3
			}
		}
	}
	my.qM(m.list)
}
#' @export

my.qM.diag.real <- function(r){
	n <- length(r)
	m.list <- list()
	for(i in 1:4){
		m.list[[i]] <- Matrix(0,n,n)
		if(i==1){
			diag(m.list[[i]]) <- r
		}
	}
	my.qM(m.list)
}

#' @export

my.D_rho <- function(mesh,xyz,area.norm=my.area.norm(mesh,xyz),face.rho){
	D <- my.Dirac(mesh,xyz,area.norm)
	P <- my.qM.diag.real(face.rho)
	B <- my.B(mesh)
	D-P%*%B
}
#' @export

my.dec.adjoint <- function(area.norm,X,noMv=TRUE){
	Mf <- my.qM.diag.real(area.norm$face.area)
	#solve(Mv) %*% X %*% Mf # 本来はこれだが、すべての計算にsolve(Mv)をかぶせるので省略できる

	if(!noMv){
		Mv <- my.qM.diag.real(area.norm$vertex.area)
	}
	if(noMv){
		return(t(X) %*% Mf)
	}else{
		return(solve(Mv) %*% t(X) %*% Mf)
	}
}
#' @export

my.V <- function(area.norm){
	my.qM.diag.real(area.norm$vertex.area^(-0.5))
}
#' @export

my.AVsq <- function(A,area.norm,noMv=TRUE){
	V <- my.V(area.norm)
	AV <- A %*% V
	my.dec.adjoint(area.norm,AV,noMv) %*% AV
}
#' @export

my.inv.pow.2 <- function(A,n.iter=3,b=rep(1,ncol(A)),log=FALSE){
  x <- b
  x <- x/sqrt(sum(x^2))
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

my.smallest.lambda <- function(A,area.norm){
	AVsq <- my.AVsq(A,area.norm)
	V <- my.V(area.norm)
	tmp.lambda <- my.inv.pow.2(AVsq)[[1]]
	V %*% tmp.lambda
}
#' @export

my.update.edge <- function(mesh,area.norm,lambda){
	e.vector.q <- Hi*area.norm$edge.vector[,1] + Hj*area.norm$edge.vector[,2] + Hk*area.norm$edge.vector[,3]
	lambda.m <- matrix(lambda,byrow=TRUE,ncol=4)
	lambda.q <- lambda.m[,1] + lambda.m[,2]*Hi + lambda.m[,3]*Hj + lambda.m[,4]*Hk
	lambda.q.bar <- lambda.m[,1] -(lambda.m[,2]*Hi + lambda.m[,3]*Hj + lambda.m[,4]*Hk)
	e.he <- mesh$E$he
	e.he.nxt <- mesh$HE$nxt[e.he]
	v.st <- mesh$HE$v[e.he]
	v.end <- mesh$HE$v[e.he.nxt]
	ret <- 1/3 * lambda.q.bar[v.st] * e.vector.q * lambda.q[v.st] + 1/6 * lambda.q.bar[v.end] * e.vector.q * lambda.q[v.st] + 1/6 * lambda.q.bar[v.end] * e.vector.q * lambda.q[v.st] + 1/3 * lambda.q.bar[v.end] * e.vector.q * lambda.q[v.end]
	t(as.matrix(ret)[2:4,])
}
#' @export

my.update.f <- function(new.e,d0,star0,star1){
	edge.q <- Hi * new.e[,1] + Hj * new.e[,2] + Hk * new.e[,3]
	L <- my.Laplacian2(d0,star0,star1)
	L.list <- list()
	L.list[[1]] <- L
	L.list[[2]] <- L.list[[3]] <- L.list[[4]] <- Matrix(0,length(L[,1]),length(L[1,]))
	L.q <- my.qM(L.list)
	
	Div <- my.div(d0,star1)
	Div.list <- list()
	Div.list[[1]] <- Div
	Div.list[[2]] <- Div.list[[3]] <- Div.list[[4]] <- Matrix(0,length(Div[,1]),length(Div[1,]))
	Div.q <- my.qM(Div.list)
	
	solve(L.q) %*% Div.q %*% c(as.matrix(edge.q))
}
#' @export

my.deform.k.serial2 <- function(mesh,xyz,area.norm,vertices,faces.v,k,n.iter=10){
  v.list <- rho.cot.list <- list()
  v.list[[1]] <- vertices
  rho.cot.list[[1]] <- my.curvature.cot(v.list[[1]],faces.v)
  for(i in 1:n.iter){
    
    tmp.rho.f <- rho.cot.list[[i]][[3]] * Mod(Im(rho.cot.list[[i]][[2]]))
D_rho <- my.D_rho(mesh,xyz,area.norm,face.rho=k*tmp.rho.f)

V <- my.V(area.norm)

AVsq <- my.AVsq(D_rho,area.norm)

lambda <- my.smallest.lambda(D_rho,area.norm)

new.e <- my.update.edge(mesh,area.norm,lambda)

new.f <- my.update.f(new.e,dec$d0,dec$star0,dec$star1)

    #tmp.out <- my.conformal.rho(v.list[[i]],faces.v,k*tmp.rho.f,face=TRUE)
    v.list[[i+1]] <- tmp.out$xyz.new.q
    rho.cot.list[[i+1]] <- my.curvature.cot(v.list[[i+1]],faces.v)
  }
  return(list(v.list=v.list,rho.cot.list=rho.cot.list,k=k))
}

#' @export

my.halfedge <- function(xyz,faces.v){
	n.v <- max(c(faces.v))
	n.f <- length(faces.v[1,])
	n.e <- n.f*3/2
	n.he <- n.e * 2
	
	HE <- list(id=1:n.he,v=rep(NA,n.he),f=rep(NA,n.he),e=rep(NA,n.he),flip=rep(NA,n.he),nxt=rep(NA,n.he))
	
	F <- list(id=1:n.f,he=rep(NA,n.f),v=matrix(NA,n.f,3),norm=matrix(NA,n.f,3),area=rep(NA,n.f))
	V <- list(id=1:n.v,he=rep(NA,n.v),f=list(),xyz=xyz,norm=matrix(NA,n.v,3),area=rep(NA,n.v))
	E <- list(id=1:n.e,he=rep(NA,n.e),cot=rep(NA,n.e))
	
	for(i in 1:n.f){
		cnt1 <- (i-1) * 3 + 1
		cnt2 <- (i-1) * 3 + 2
		cnt3 <- (i-1) * 3 + 3
		HE$v[cnt1] <- faces.v[1,i]
		HE$v[cnt2] <- faces.v[2,i]
		HE$v[cnt3] <- faces.v[3,i]
		HE$f[c(cnt1,cnt2,cnt3)] <- i
		HE$nxt[cnt1] <- cnt2
		HE$nxt[cnt2] <- cnt3
		HE$nxt[cnt3] <- cnt1
		F$he[i] <- cnt1
		if(is.na(V$he[faces.v[1,i]])){
			V$he[faces.v[1,i]] <- cnt1
		}
		if(is.na(V$he[faces.v[2,i]])){
			V$he[faces.v[2,i]] <- cnt2
		}
		if(is.na(V$he[faces.v[3,i]])){
			V$he[faces.v[3,i]] <- cnt3
		}
	}
	
	v.st <- v.end <- rep(NA,n.he)
	for(i in 1:n.he){
		v.st[i] <- HE$v[i]
		v.end[i] <- HE$v[HE$nxt[i]]
	}
	v.st.end <- t(apply(cbind(v.st,v.end),1,sort))
	v.st.end.value <- v.st.end[,1]-1 + v.st.end[,2]*(n.v-1)
	ord <- order(v.st.end.value)
	ord2 <- matrix(ord,byrow=TRUE,ncol=2)
	HE$flip[ord2[,1]] <- ord2[,2]
	HE$flip[ord2[,2]] <- ord2[,1]
	HE$e[ord2[,1]] <- 1:n.e
	HE$e[ord2[,2]] <- 1:n.e
	E$he <- ord2[,1]
	
	V$f <- my.enumerate.v.f2(V,HE)
	F$v <- my.enumerate.f.v2(F,HE)
	
	for(i in 1:n.f){
		tmp <- my.f.area.norm(V$xyz[F$v[i,],])
		F$area[i] <- tmp$area
		F$norm[i,] <- tmp$norm
	}
	for(i in 1:n.v){
		V$area[i] <- sum(F$area[V$f[[i]]])/3
		tmp <- apply(F$norm[V$f[[i]],],2,sum)
		V$norm[i,] <- tmp/sqrt(sum(tmp^2))
	}
	
	for(i in 1:n.he){
		p0 <- HE$v[HE$nxt[HE$nxt[i]]]
		p1 <- HE$v[i]
		p2 <- HE$v[HE$nxt[i]]
		v1 <- V$xyz[p1,]-V$xyz[p0,]
		v2 <- V$xyz[p2,]-V$xyz[p0,]
		tmp <- c(v1[2]*v2[3]-v1[3]*v2[2],v1[3]*v2[1]-v1[1]*v2[3],v1[1]*v2[2]-v1[2]*v2[1])
		HE$cot[i] <- sum(v1*v2)/sqrt(sum(tmp^2))
		
	}
	star0 <- Diagonal(x=V$area)
	e.cot <- rep(NA,n.e)
	d0 <- Matrix(0,n.e,n.v)
	for(i in 1:n.e){
		he <- E$he[i]
		flip.he <- HE$flip[he]
		e.cot[i] <- (HE$cot[he] + HE$cot[flip.he])/2
		d0[i,HE$v[he]] <- -1
		d0[i,HE$v[flip.he]] <- 1
	}
	star1 <- Diagonal(x=e.cot)
	star2 <- Diagonal(x=F$area)

	Laplacian <- solve(star0) %*% t(d0) %*% star1 %*% d0

	return(list(HE=HE,V=V,E=E,F=F,star0=star0,star1=star1,d0=d0,Laplacian=Laplacian))
	
}
#' @export

my.f.area.norm <- function(X){
	p0 <- X[1,]
	p1 <- X[2,]
	p2 <- X[3,]
	v1 <- p1-p0
	v2 <- p2-p0
	tmp <- c(v1[2]*v2[3]-v1[3]*v2[2],v1[3]*v2[1]-v1[1]*v2[3],v1[1]*v2[2]-v1[2]*v2[1])
	L <- sqrt(sum(tmp^2))
	return(list(area=L/2,norm=tmp/L))
}
#' @export
my.enumerate.v.he <- function(A,vid){
	ret <- list()
	for(j in 1:length(vid)){
		i <- vid[j]
		this.he <- A$V$he[i]
		tmp.he <- A$HE$flip[this.he]
		next.he <- A$HE$nxt[tmp.he]
		ret[[j]] <- c(next.he)
		while(this.he!=next.he){
			tmp.he <- A$HE$flip[next.he]
			next.he <- A$HE$nxt[tmp.he]
			ret[[j]] <- c(ret[[j]],next.he)
		}
	}
	ret
}

#' @export

# 頂点を含む面を列挙する
my.enumerate.v.f <- function(A,vid){
	ret <- list()
	for(j in 1:length(vid)){
		i <- vid[j]
		this.he <- A$V$he[i]
		tmp.he <- A$HE$flip[this.he]
		next.he <- A$HE$nxt[tmp.he]
		ret[[j]] <- c(A$HE$f[next.he])
		while(this.he!=next.he){
			tmp.he <- A$HE$flip[next.he]
			next.he <- A$HE$nxt[tmp.he]
			ret[[j]] <- c(ret[[j]],A$HE$f[next.he])
		}
	}
	ret
}
#' @export

my.enumerate.v.f2 <- function(V,HE){
	ret <- list()
	vid <- V$id
	for(j in 1:length(vid)){
		i <- vid[j]
		this.he <- V$he[i]
		tmp.he <- HE$flip[this.he]
		next.he <- HE$nxt[tmp.he]
		ret[[j]] <- c(HE$f[next.he])
		while(this.he!=next.he){
			tmp.he <- HE$flip[next.he]
			next.he <- HE$nxt[tmp.he]
			ret[[j]] <- c(ret[[j]],HE$f[next.he])
		}
	}
	ret
}
#' @export

# 面の３頂点を列挙する
my.enumerate.f.v <- function(A,fid){
	ret <- list()
	for(j in 1:length(fid)){
		i <- fid[j]
		this.he <- A$F$he[i]
		next.he <- A$HE$nxt[this.he]
		ret[[j]] <- c(A$HE$v[next.he])
		while(this.he!=next.he){
			next.he <- A$HE$nxt[next.he]
			ret[[j]] <- c(ret[[j]],A$HE$v[next.he])
		}
	}
	ret
}
#' @export

my.enumerate.f.v2 <- function(F,HE){
	fid <- F$id
	ret <- matrix(NA,length(fid),3)
	
	for(j in 1:length(fid)){
		i <- fid[j]
		this.he <- F$he[i]
		next.he <- HE$nxt[this.he]
		tmp <- c(HE$v[next.he])
		while(this.he!=next.he){
			next.he <- HE$nxt[next.he]
			tmp <- c(tmp,HE$v[next.he])
		}
		ret[j,] <- tmp
	}
	ret
}
#' @export

my.face.norm <- function(A,fid){
	ret <- matrix(NA,length(fid),3)
	for(j in 1:length(fid)){
		i <- fid[j]
		
	}
}

#' @export
my.back.euler <- function(xyz,Laplacian,step=0.1){
	n.v <- length(vertices)
	I <- Diagonal(n.v)
	new.x <- solve(I-step*Laplacian, xyz[,1])
	new.y <- solve(I-step*Laplacian, xyz[,2])
	new.z <- solve(I-step*Laplacian, xyz[,3])
	new.xyz <- cbind(new.x,new.y,new.z)
	new.xyz
}
