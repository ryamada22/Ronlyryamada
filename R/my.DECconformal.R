#' Utilities for DEC conformal
#'
#' These functions are for DEC conformal
#' @export


#' 頂点のrho値から面のrho値を算出する
rho.fromVtoTri <- function(rho.v,faces.v){
  tmp.rho <- matrix(rho.v[faces.v],nrow=3)
  apply(tmp.rho,2,mean)
}
#' Utility 関数
#' 特に疎行列において、面単位で不定回数の値加算をするためのユーティリティ関数
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
#' E行列作成関数
#' ただし、四元数行列なので、虚実４成分それぞれをリストの要素とし
#' 個々のリストも行列ではなく、ベクトル化した形で返す
#' 疎行列のリスト
#' 三角形の面積を考慮して、面のrho値に応じて、rhoが作る曲面上の
#' ベクトル場を作り
#' その定常解の一つとして固有ベクトルを取る
#' @export
my.make.E.v <- function(vertices,faces.v,rho){
  # 三角形の面積
	edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
	edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	tmp <- edge1 * edge2
  A <- Mod(Im(tmp))/2
	#A <- abs(i(tmp)+j(tmp)+k(tmp))/2
	# 三角形ごとの計算用係数
	coef.a <- -1/(4*A)
	coef.b <- rho/6
	coef.c <- A*rho^2/9
	
	# Rでは四元数を要素とする行列がないので、re,i,j,kごとに正方行列を作ることにする
	E.re <- E.i <- E.j <- E.k <- sparseVector(c(0),i=c(1),length=length(vertices)^2)
  
	e.q <- list()
	e.q[[1]] <- vertices[faces.v[2,]]-vertices[faces.v[3,]]
	e.q[[2]] <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
	e.q[[3]] <- vertices[faces.v[1,]]-vertices[faces.v[2,]]
  # すべての頂点ペアについて、すべての三角形ごとに処理をする
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
#' 四元数行列を表すベクトルリストを引数として
#' その実数版行列を返す
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
#' Ax = lambda b の固有ベクトルを逆冪乗法で近似する
#' 初期ベクトル bは(1,1,...)がよいことが知られているので
#' それをデフォルト値としてある
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
#' ラプラシアン行列の算出
#' すべての点について、その点を含む三角形について、三角形の内角で定まる値を加算する
#' @export
my.make.L <- function(vertices,faces.v){
  n.v <- length(vertices)
	L <- sparseVector(c(0),i=c(1),length=n.v^2)
	for(i in 1:3){
		v.ord <- ((1:3)+i+1) %% 3 + 1
		k1 <- faces.v[v.ord[1],]
		k2 <- faces.v[v.ord[2],]
		k3 <- faces.v[v.ord[3],]
		# 頂点四元数
		v1 <- vertices[k1]
		v2 <- vertices[k2]
		v3 <- vertices[k3]
		
		# edge 四元数
		u1 <- v2-v1
		u2 <- v3-v1
		# edge 四元数(純虚四元数)の積は実部がドット積、虚部がクロス積ベクトル
		u12 <- u1 * u2
		cotAlpha <- (-Re(u12))/Mod(Im(u12))
		# このcotAlphaを行列Lの[k2,k2],[k3,k3],[k2,k3],[k3,k2]に加算する
		# 疎ベクトルで格納する
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
#' ラプラシアン行列Lは実行列だが、その四元数版の実行列化したもの(行・列がそれぞれ４倍)を作る
#' @export
my.make.quatList <- function(L){
  L.q <- list()
  L.q[[1]] <- L
  for(i in 2:4){
    L.q[[i]] <- L*0
  }
  L.q
}
#' omegaは回転を定めるlambdaを介して算出する
#' @export
my.make.omega <- function(vertices,faces.v,lambda){
  n.v <- length(vertices)
	omega <- rep(0*Hi,n.v)
	for(i in 1:3){
		v.ord <- ((1:3)+i+1) %% 3 + 1
		k1 <- faces.v[v.ord[1],]
		# 対向辺の向きは頂点IDの大小順にそろえる
		k23 <- rbind(faces.v[v.ord[2],],faces.v[v.ord[3],])
		k23 <- apply(k23,2,sort)
		k2 <- k23[2,]
		k3 <- k23[1,]
		# 頂点四元数
		v1 <- vertices[k1]
		v2 <- vertices[k2]
		v3 <- vertices[k3]
		
		edge <- v3-v2
		
		# lambdaの四元数化、とその共役四元数化
		lambda.mat <- matrix(lambda,nrow=4)
		lambda.q <- as.quaternion(lambda.mat)
		lambda.mat.2 <- lambda.mat
		lambda.mat.2[2:4,] <- -lambda.mat.2[2:4,]
		lambda.q. <- as.quaternion(lambda.mat.2)
		lambda1 <- lambda.q[k2]
		lambda1. <- lambda.q.[k2]
		lambda2 <- lambda.q[k3]
		lambda2. <- lambda.q.[k3]
		# エッジに回転処理をしながら、加算量を算出
		val <- 1/3 * lambda1. * edge * lambda1 + 1/6 * lambda1. * edge * lambda2 + 1/6 * lambda2. * edge * lambda1 + 1/3 * lambda2. * edge * lambda2

		# edge 四元数
		u1 <- v2-v1
		u2 <- v3-v1
		# edge 四元数(純虚四元数)の積は実部がドット積、虚部がクロス積ベクトル
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
#' 局所平均曲率を計算する関数
#' 入力は頂点の四元数座標ベクトルと３行で表された三角形頂点ID行列
#' 返り値はnorm.v,norm.face=ret.face,dirの３つで
#' norm.vは各頂点の法線ベクトルに相当する四元数の行列で、その虚部の絶対値が平均曲率の大きさ
#' norm.faceは各三角形のそれで、虚部の絶対値が平均曲率の大きさ
#' dirは各三角形の平均曲率が曲面の外向きか内向きか(今、扱っている曲面はすべて球面様の閉曲面なので、外向き・内向きはその意味での内外)を示すc(-1,0,1)
#' @export
my.curvature.cot <- function(vertices,faces.v){
  n.v <- length(vertices)
  ret <- rep(0*Hk,n.v)
  # 三角形の面積
  edge1 <- vertices[faces.v[2,]]-vertices[faces.v[1,]]
  edge2 <- vertices[faces.v[3,]]-vertices[faces.v[1,]]
  tmp <- edge1 * edge2
  tmp2 <- Mod(Im(tmp))
  #tmp2 <- i(tmp)+j(tmp)+k(tmp)
  #inv.f <- which(tmp2<0)
  #faces.v[2:3,inv.f] <- rbind(faces.v[3,inv.f],faces.v[2,inv.f])
 	A.f <- (tmp2)/2
  # 頂点周囲の三角形面積の和
  val <- rep(A.f,each=3)
  addr <- c(faces.v)
  tmp.out <- my.vector.access(val,addr)
  A.v <- tmp.out[[2]][,1] # 頂点周囲面積和
  
  for(i in 1:3){
		v.ord <- ((1:3)+i+1) %% 3 + 1
		k1 <- faces.v[v.ord[1],]
		k2 <- faces.v[v.ord[2],]
		k3 <- faces.v[v.ord[3],]
		# 頂点四元数
		v1 <- vertices[k1]
		v2 <- vertices[k2]
		v3 <- vertices[k3]
		
		# edge 四元数
		u1 <- v2-v1
		u2 <- v3-v1
    
    # 対向辺
    u3 <- v3-v2
    #u3.len <- Mod(Im(u3))
		# edge 四元数(純虚四元数)の積は実部がドット積、虚部がクロス積ベクトル
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

#' ３列の実行列(３次元座標の行列)xyzと、３行の頂点ID行列faces.vとを用いて、三角形メッシュを塗り描く関数
#' ただし、この関数は、元がsp.triは、元の三角形メッシュがmy.sphere.tri.mesh()によって作られていることを前提としている。そのような場合に、sp.triという引数によって、隣接緯線間ごとに色を塗り分けることで、元の正球からの変形の様子を色表現することができる
#' また、各三角形の色はrho.fなる実ベクトルによってコントラストをつけることができる
#' col1は、元の正球からの変化を縞模様化するための長さ２の整数ベクトルである。その縞模様を使いたくないとき(メッシュがmy.sphere.tri.mesh()由来ではないときを含む)は、この長さ２のベクトルの２つの値を同じにすればよい)
#' @export
plot.sp.conformal <- function(xyz,faces.v,sp.tri,rho.f,col1=c(4,5)){
  plot3d(xyz,xlab="x",ylab="y",zlab="z")
  mesh.tri <- tmesh3d(t(xyz),faces.v,homogeneous=FALSE)
  
  # 縞のための値ベクトル
  col. <- rep(col1,length(sp.tri))[1:length(sp.tri)]
  col <- rep(col.,sp.tri*3)
  # rhoを反映した値ベクトル
  rho.f <- rep(rho.f,each=3)
  #rho.f2 <- sign(rho.f)
  rho.f <- (rho.f-min(rho.f))/(max(rho.f)-min(rho.f))
  
  col2 <- rgb(1-rho.f,1,col/6)
  #col3 <- gray((rho.f2+1)*0.5)
  shade3d(mesh.tri,col=col2)  
}