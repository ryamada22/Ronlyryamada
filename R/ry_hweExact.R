#' Hardy-Weinberg Exact Test
#'
#' This function returns Exact HWtest p and probability of all possible counts.
#' @param g is a integer vector with length 3.
#' @keywords exact test, Hardy-Weinberg equilibrium
#' @export
#' @examples
#' xx<-ry_hweExact(c(813,182,5))
#' xx$p.value
#' sum(xx$prob) # should be 1

ry_hweExact<-function(g=c(813,182,5)){
 n<-sum(g) # total indvidual number
 nA<-2*g[1]+g[2] # No. A-alleles
 na<-2*g[3]+g[2] # No. a-alleles
 evod<-g[2]%%2
 maxAa<-min(nA,na)-evod
 Aa<-seq(from=evod,to=maxAa,by=2)
 AA<-(nA-Aa)/2 
 aa<-(na-Aa)/2 
 obs<-(g[2]-evod)/2+1 
 prob<-rep(0,length(Aa)) 
 prob<-exp(n*lgamma(2+1)+lgamma(nA+1)+lgamma(na+1)-lgamma(2*n+1)-(AA*lgamma(2+1)+Aa*lgamma(1+1)+aa*lgamma(2+1))+lgamma(n+1)-(lgamma(AA+1)+lgamma(Aa+1)+lgamma(aa+1)))
 p.value<-sum(prob[prob<=prob[obs]])
 
 list(Aa=Aa,prob=prob,obsprob=prob[obs],p.value=p.value)
}

