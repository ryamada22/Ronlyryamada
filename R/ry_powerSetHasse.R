#' Draw a Hasse Diagram
#'
#' This function draws a Hasse diagram..
#' @param k is an integer..
#' @keywords Hasse powerSet
#' @export
#' @examples
#' ry_powerSetHasse(4)

ry_powerSetHasse<-function(k=4){
 lev<-k+1
 numelem<-rep(0,lev)
 a<-0:k
 numelem<-choose(k,a)

 numnode<-2^k
 nodesX<-rep(0,numnode)
 nodesY<-rep(0,numnode)
 vecs<-matrix(rep(0,numnode*k),ncol=k)
 keta<-2^(0:(k-1))

 counter<-rep(0,lev)
 startval<--(numelem-1)/2
 for(i in 2:numnode){
  vecs[i,]<-vecs[i-1,]
  vecs[i,1]<-vecs[i,1]+1
  for(j in 1:(k-1)){
   if(vecs[i,j]==2){
    vecs[i,j]<-0
    vecs[i,j+1]<-vecs[i,j+1]+1
   }else{
    j<-k
   }
  }
  tmpsum<-sum(vecs[i,])
  nodesY[i]<-tmpsum
  nodesX[i]<-counter[tmpsum]+startval[tmpsum+1]
  counter[tmpsum]<-counter[tmpsum]+1
 }

 values<-vecs%*%keta

 edges<-matrix(rep(FALSE,numnode^2),ncol=numnode)
 for(i in 1:(numnode-1)){
  for(j in (i+1):numnode){
   tmpdist<-sum(abs(vecs[i,]-vecs[j,]))
   if(tmpdist==1){
    edges[i,j]<-TRUE
   }
  }
 }
 Edges<-which(edges,arr.ind=TRUE)
 plot(nodesX,nodesY,pch=15)
 segments(nodesX[Edges[,1]],nodesY[Edges[,1]],nodesX[Edges[,2]],nodesY[Edges[,2]])
}
