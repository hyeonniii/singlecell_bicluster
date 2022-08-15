library(Rcpp)
library(gtools)
library(Matrix)
library(dplyr)
library(parallel)
library(doSNOW)
library(aricode)

numCores<-4
cl<-makeCluster(numCores)
registerDoSNOW(cl)
library(foreach)

setwd("C:\\sc_bicluster")
sourceCpp("scfn_revised.cpp")

set.seed(123)

n<-40
p<-80
K<-3
R<-4
truthCs<-sample(1:K,n,rep=TRUE)
truthDs<-sample(1:R,p,rep=TRUE)
mus<-runif(K*R,-2,2)

mus<-matrix(c(mus),nrow=K,ncol=R,byrow=F)

x<-matrix(rnorm(n*p,mean=0,sd=1),nrow=n,ncol=p)

musmatrix<-matrix(NA,nrow=n,ncol=p)
for(i in 1:max(truthCs)){
  for(j in 1:max(truthDs)){
    x[truthCs==i,truthDs==j]<-x[truthCs==i,truthDs==j]+mus[i,j]
    musmatrix[truthCs==i,truthDs==j]<-mus[i,j]
  }
}

gs<-list(lambda=c(5,6,6.5,7),gamma=c(0.5,1)) %>% expand.grid()



idx<-sort(sample(1:(n*p),n*p*0.2))
vecx<-c(x)
vecx[idx]<-0
x<-matrix(vecx,nrow=n,ncol=p)

rs2<-foreach(i=1:nrow(gs),.packages='Rcpp',.noexport="sc_bicluster")%dopar%{
  sourceCpp("C:\\sc_bicluster\\scfn_revised.cpp")
  result<-sc_bicluster(as.numeric(gs[i,1]),as.numeric(gs[i,2]),0.5,x,x,0.01,50)
  list(X=result$X,iter_error=result$iter_error)
}

km.Cs<-kmeans(x,K,nstart=20)$cluster
km.Ds<-kmeans(t(x),R,nstart=20)$cluster
km.mus<-matrix(NA,nrow=n,ncol=p)


for(i in 1:n){
  for(j in 1:p){
    km.mus[i,j]<-mean(x[km.Cs==km.Cs[i],km.Ds==km.Ds[j]])
  }
}

rs<-foreach(i=1:nrow(gs),.packages='Rcpp',.noexport="sc_bicluster")%dopar%{
  sourceCpp("C:\\sc_bicluster\\scfn_revised.cpp")
  result<-sc_bicluster(as.numeric(gs[i,1]),as.numeric(gs[i,2]),0.5,x,km.mus,0.01,50)
  list(X=result$X,iter_error=result$iter_error)
}


gs<-list(lambda=c(3,3.5,4,4.5),gamma=c(0.5,1)) %>% expand.grid()

rs3<-foreach(i=1:nrow(gs),.packages='Rcpp',.noexport="sc_bicluster")%dopar%{
  sourceCpp("C:\\sc_bicluster\\scfn_revised.cpp")
  result<-sc_bicluster(as.numeric(gs[i,1]),as.numeric(gs[i,2]),0.5,x,km.mus,0.01,50)
  list(X=result$X,iter_error=result$iter_error)
}

gs<-list(lambda=c(3,3.5,4,4.5,5,6),gamma=c(0.1,0.5,1)) %>% expand.grid()

##soarse하지 않은 세팅..
rs4<-foreach(i=1:nrow(gs),.packages='Rcpp',.noexport="sc_bicluster")%dopar%{
  sourceCpp("C:\\sc_bicluster\\scfn_revised.cpp")
  result<-sc_bicluster(as.numeric(gs[i,1]),as.numeric(gs[i,2]),0.5,x,x,0.01,50)
  list(X=result$X,iter_error=result$iter_error)
}

#다시 sparse한 세팅

idx<-sort(sample(1:(n*p),n*p*0.2))
vecx<-c(x)
vecx[idx]<-0
x<-matrix(vecx,nrow=n,ncol=p)

gs<-list(lambda=c(1,1.5,2,2.5),gamma=c(0)) %>% expand.grid()

rs<-foreach(i=1:nrow(gs),.packages='Rcpp',.noexport="sc_bicluster")%dopar%{
  sourceCpp("C:\\sc_bicluster\\scfn_revised.cpp")
  result<-sc_bicluster(as.numeric(gs[i,1]),as.numeric(gs[i,2]),0.5,x,x,0.01,50)
  list(X=result$X,iter_error=result$iter_error)
}

stopCluster(cl)

#############

library(aricode)
nmifn<-function(x,k,reallabel){
  distmat<-dist(x)
  h1<-hclust(distmat)
  labels<-cutree(h1,k=k)
  NMI(reallabel,labels)
}

for (idx in 1:nrow(gs)){
  print(nmifn(rs5[[idx]][["X"]],K,truthCs))
}

nmifn(x,K,truthCs)
NMI(km.Cs,truthCs)
