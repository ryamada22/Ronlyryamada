#' Euler Triangulation on Sphere
#'
#' A set of functions to randomly generate euler triangulation on Sphere.
#' @examples
#' library(igraph)
#' n.step <- 20
#' euler.out <- my.EulTriSph(n.step)
#' el <- my.dir.graph(euler.out[[1]],euler.out[[2]])
#' plot(graph.edgelist(el))
#' root.e <- el[1,]
#' euler.tree.out <- my.euler.tree(euler.out[[1]],euler.out[[2]],root.e)
#' g.tree <- graph.edgelist(euler.tree.out[[1]])
#' is.connected(g.tree)
#' g.tree
#' sum(euler.out[[2]])
#' nrow(euler.out[[1]])
#' sh.dist <- shortest.paths(graph.edgelist(euler.tree.out[[2]]),root.e[1],mode="out")
#' plot(graph.edgelist(euler.tree.out[[1]]),vertex.label=sh.dist)
#' @import igraph
#' @export
my.EulTriSph <- function(n.step){
  # エッジリスト
  el <- matrix(c(1,2,1,3,1,4,1,5,2,3,3,4,4,5,5,2,2,6,3,6,4,6,5,6),byrow=TRUE,ncol=2)
  # 頂点座標
  x <- matrix(c(0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1),byrow=TRUE,ncol=3)
  x <- x/sqrt(apply(x^2,1,sum))
  # 隣接頂点、順序考慮のリスト
  eofv <- list()
  eofv[[1]] <- c(2,3,4,5)
  eofv[[2]] <- c(1,5,6,3)
  eofv[[3]] <- c(1,2,6,4)
  eofv[[4]] <- c(1,3,6,5)
  eofv[[5]] <- c(1,4,6,2)
  eofv[[6]] <- c(5,4,3,2)
  
  # 三角形構成頂点の行列
  # 三角形のノードのIDでソートしておく
  tris <- matrix(c(1,2,3,1,3,4,1,4,5,1,5,2,6,2,3,6,3,4,6,4,5,6,5,2),byrow=TRUE,ncol=3)
  tris <- t(apply(tris,1,sort))
  
  # 三角形の向きをそろえておく
  tris.ord <- matrix(c(1,2,3,1,4,3,1,4,5,1,2,5,6,2,3,6,4,3,6,4,5,6,2,5),byrow=TRUE,ncol=3)
  
  # 三角形の色
  cols <- c(FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE)

  # 処理回数を指定して実行
  for(i in 1:n.step){
    # p/q処理の選択と対象頂点を選択
    proc <- my.process.choice(eofv)

  	if(proc$type==1){# p処理
  		# rule p
  		out <- rule.p(eofv,proc$v,x,tris,cols,tris.ord)
  		eofv <- out[[1]]
  		x <- out[[2]]
  		tris <- out[[3]]
  		cols <- out[[4]]
  		tris.ord <- out[[5]]
  	}else{# q処理
  		# rule q
  		out <- rule.q(eofv,proc$v,x,tris,cols,tris.ord)
  		eofv <- out[[1]]
  		x <- out[[2]]
  		tris <- out[[3]]
  		cols <- out[[4]]
  		tris.ord <- out[[5]]
  	}
  }
  return(list(tris.ord=tris.ord,cols=cols,tris=tris,eofv=eofv,x=x))
}

#' @export
rule.p <- function(eofv,v,x,tris,cols,tris.ord){
  # 処理対象頂点vのどの部分にひし形を挿入するかをランダムに決める
  st.v <- sample(1:length(eofv[[v]]),1)
  # 扱う隣接頂点のリスト
	nb <- c(eofv[[v]],eofv[[v]])[st.v:(st.v+2)]
	# 新規頂点のID
	newv <- 1:2 + length(eofv)
	
  # 隣接頂点ベクトルの更新
	eofv[[nb[1]]] <- my.replace.nodes(eofv[[nb[1]]],v,c(newv[2],newv[1],v))
	eofv[[nb[2]]] <- my.replace.nodes(eofv[[nb[2]]],v,c(newv[2]))
	eofv[[nb[3]]] <- my.replace.nodes(eofv[[nb[3]]],v,c(v,newv[1],newv[2]))
	eofv[[v]] <- my.replace.nodes(eofv[[v]],nb[2],c(newv[1]))
	
	eofv[[newv[1]]] <- c(v,nb[1],newv[2],nb[3])
	eofv[[newv[2]]] <- c(newv[1],nb)
	
  # 三角形の色、辺の向きを定める
  # 定めるにあたり、処理対称部分の色パターン、辺の向きパターンを確認
	tmp.tri <- sort(c(v,nb[1],nb[2]))
	tmp.tri.id <- which(tris[,1]==tmp.tri[1] & tris[,2] ==tmp.tri[2] & tris[,3]==tmp.tri[3])
	tmp.col <- cols[tmp.tri.id]
  
	add.tris <- matrix(c(v,nb[1],newv[1],v,nb[3],newv[1],newv[1],newv[2],nb[1],newv[1],newv[2],nb[3],newv[2],nb[1],nb[2],newv[2],nb[3],nb[2]),byrow=TRUE,ncol=3)
	add.tris.ord <- add.tris

  # 処理対象領域の状況に合わせて更新情報を改変
	tmp.order.tri <- tris.ord[tmp.tri.id,]
  tmp.order.tri2 <- c(tmp.order.tri,tmp.order.tri[1])
  loc1 <- which(tmp.order.tri==v)
  if(tmp.order.tri2[loc1+1]==nb[1]){
    new.tris.ord <- rbind(tris.ord,add.tris.ord)
  }else{
    new.tris.ord <- rbind(tris.ord,cbind(add.tris.ord[,3],add.tris.ord[,2],add.tris.ord[,1]))
  }

  # 三角形構成頂点ID(値ソート)、色ベクトルの更新
	add.tris <- t(apply(add.tris,1,sort))
	add.cols <- c(TRUE,FALSE,FALSE,TRUE,TRUE,FALSE)
	if(!tmp.col){
		add.cols <- !add.cols
	}
	newtris <- rbind(tris,add.tris)
	newcols <- c(cols,add.cols)
	
  # 処理によって除去される三角形についての処理
	rmtris1 <- sort(c(v,nb[1],nb[2]))
	rmtris2 <- sort(c(v,nb[2],nb[3]))
	
	rmtris1id <- which(newtris[,1]==rmtris1[1] & newtris[,2]==rmtris1[2] & newtris[,3]==rmtris1[3])
	rmtris2id <- which(newtris[,1]==rmtris2[1] & newtris[,2]==rmtris2[2] & newtris[,3]==rmtris2[3])
	
	newtris <- newtris[-c(rmtris1id,rmtris2id),]
	newcols <- newcols[-c(rmtris1id,rmtris2id)]
	
	new.tris.ord <- new.tris.ord[-c(rmtris1id,rmtris2id),]
	
  # 座標の更新
	newx <- rbind(x,matrix(0,2,3))
	id1 <- which(eofv[[v]]==nb[1])
	id2 <- which(eofv[[v]]==nb[2])
	id3 <- which(eofv[[v]]==nb[3])
	rest <- eofv[[v]][-(c(id1,id2,id3))]
	#newx[v,] <- (apply(matrix(x[c(id1,id2,id3),],ncol=3),2,mean) + x[v,])/2
	newx[newv[1],] <- (2*newx[v,] + newx[nb[2],])/3

	newx[newv[2],] <- (newx[v,] + 2 * newx[nb[2],])/3
	newx <- newx/sqrt(apply(newx^2,1,sum))
	return(list(eofv=eofv,x=newx,tris=newtris,cols=newcols,tris.ord=new.tris.ord))
}

#' @export
rule.q <- function(eofv,v,x,tris,cols,tris.ord){
	loop <- TRUE
  # 周囲ノードでnb[2],nb[5]として取り出された２頂点間にすでに辺があるときは
  # 同一の処理ができないのでそうならない場合を取り出す
	while(loop){
		st.v <- sample(1:length(eofv[[v]]),1)
		nb <- c(eofv[[v]],eofv[[v]])[st.v:(st.v+4)]
		if(length(which(eofv[[nb[2]]]==nb[5])) == 0){
			loop <- FALSE
		}
	}

	newv <- 1 + length(eofv)
	
	eofv[[nb[2]]] <- my.replace.nodes(eofv[[nb[2]]],v,c(newv,nb[5],v))
	eofv[[nb[3]]] <- my.replace.nodes(eofv[[nb[3]]],v,c(newv))
	eofv[[nb[4]]] <- my.replace.nodes(eofv[[nb[4]]],v,c(newv))
	eofv[[nb[5]]] <- my.replace.nodes(eofv[[nb[5]]],v,c(v,nb[2],newv))
	
	take1 <- which(eofv[[v]]==nb[3])
	take2 <- which(eofv[[v]]==nb[4])
	eofv[[v]] <- eofv[[v]][-c(take1,take2)]
	
	eofv[[newv]] <- nb[2:5]
	
	tmp.tri <- sort(c(v,nb[2],nb[3]))
	tmp.tri.id <- which(tris[,1]==tmp.tri[1] & tris[,2] ==tmp.tri[2] & tris[,3]==tmp.tri[3])
	tmp.col <- cols[tmp.tri.id]
	
	add.tris <- matrix(c(v,nb[2],nb[5],nb[2],nb[5],newv,newv,nb[2],nb[3],newv,nb[4],nb[3],newv,nb[4],nb[5]),byrow=TRUE,ncol=3)
	add.tris.ord <- add.tris
	
  # 処理対象領域の状況に合わせて更新情報を改変
  tmp.order.tri <- tris.ord[tmp.tri.id,]
	tmp.order.tri2 <- c(tmp.order.tri,tmp.order.tri[1])
  loc1 <- which(tmp.order.tri==v)
  if(tmp.order.tri2[loc1+1]==nb[2]){
    new.tris.ord <- rbind(tris.ord,add.tris.ord)
  }else{
    new.tris.ord <- rbind(tris.ord,cbind(add.tris.ord[,3],add.tris.ord[,2],add.tris.ord[,1]))
  }

	add.tris <- t(apply(add.tris,1,sort))

	add.cols <- c(TRUE,FALSE,TRUE,FALSE,TRUE)
	if(!tmp.col){
		add.cols <- !add.cols
	}
	newtris <- rbind(tris,add.tris)
	newcols <- c(cols,add.cols)

	rmtris1 <- sort(c(v,nb[2],nb[3]))
	rmtris2 <- sort(c(v,nb[3],nb[4]))
	rmtris3 <- sort(c(v,nb[4],nb[5]))
	
	rmtris1id <- which(newtris[,1]==rmtris1[1] & newtris[,2]==rmtris1[2] & newtris[,3]==rmtris1[3])
	rmtris2id <- which(newtris[,1]==rmtris2[1] & newtris[,2]==rmtris2[2] & newtris[,3]==rmtris2[3])
	rmtris3id <- which(newtris[,1]==rmtris3[1] & newtris[,2]==rmtris3[2] & newtris[,3]==rmtris3[3])
	
	newtris <- newtris[-c(rmtris1id,rmtris2id,rmtris3id),]
	newcols <- newcols[-c(rmtris1id,rmtris2id,rmtris3id)]

	new.tris.ord <- new.tris.ord[-c(rmtris1id,rmtris2id,rmtris3id),]

	newx <- rbind(x,matrix(0,1,3))
	#rest <- eofv[[v]][-((st.v+1):(st.v+4))]
	id1 <- which(eofv[[v]]==nb[2])
	id2 <- which(eofv[[v]]==nb[3])
	id3 <- which(eofv[[v]]==nb[4])
	id4 <- which(eofv[[v]]==nb[5])
	
	rest <- eofv[[v]][-(c(id2,id3))]

	#newx[v,] <- (apply(matrix(x[rest,],ncol=3),2,mean) + x[v,])/2
	#newx[newv[1],] <- (apply(x[nb[3:4],],2,mean) + x[v,])/2
	newx[newv[1],] <- (newx[nb[3],]+newx[nb[3],])/2
	newx <- newx/sqrt(apply(newx^2,1,sum))

	return(list(eofv=eofv,x=newx,tris=newtris,cols=newcols,tris.ord=new.tris.ord))
}

#' @export
my.process.choice <- function(eofv){
  current.deg <- sapply(eofv,length)
	p.cand <- which(current.deg>=4)
	q.cand <- which(current.deg>=6)
	
	s <- sample(1:length(c(p.cand,q.cand)),1)
	type <- 0
	v <- 0
	if(s <= length(p.cand)){
		type <- 1
		v <- p.cand[s]
	}else{
		type <- 2
		v <- q.cand[s - length(p.cand)]
	}
	return(list(type=type,v=v))
}
#' @export
my.replace.nodes <- function(vec,v,u){
	# vecにあるvの位置をuに置き換える
	loc <- which(vec==v)
	pre <- c()
	post <- c()
	if(loc==1){
		post <- vec[-loc]
	}else if(loc==length(vec)){
		pre <- vec[-loc]
	}else{
		pre <- vec[1:(loc-1)]
		post <- vec[(loc+1):length(vec)]
	}
	c(pre,u,post)
}

#' @export
my.dir.graph <- function(tris,cols){
  tris.half <- tris[which(cols==0),]
  rbind(cbind(tris.half[,1],tris.half[,2]),cbind(tris.half[,2],tris.half[,3]),cbind(tris.half[,3],tris.half[,1]))
}

#' @export
my.euler.tree <- function(tris,cols,root.e){
  el <- my.dir.graph(tris,cols)
  root.st <- root.e[1]
	root.end <- root.e[2]
	rm.e <- which(el[,1]==root.st | el[,2]==root.st)

	g <- graph.edgelist(el)
  sh.dist <- shortest.paths(g,root.st,mode="out")
  
  tris.half <- tris[which(cols==0),]
  
  tris.label <- matrix(sh.dist[c(tris.half)],ncol=3)
  ord.tris.label <- t(apply(tris.label,1,order))
  
  selected.1 <- t(tris.half)[3*(0:(length(tris.half[,1])-1))+ord.tris.label[,2]]
  selected.2 <- t(tris.half)[3*(0:(length(tris.half[,1])-1))+ord.tris.label[,3]]
  
  selected.edges <- cbind(selected.1,selected.2)
  selected.edges <- rbind(root.e,selected.edges)

	list(tree.edge=selected.edges,whole.edge=el,root.e=root.e)
}


