
library(gtools)
library(dplyr)

X<-x
x<-km.mus
lambda<-0.01
gamma<-0.01
#del<-1 0.5 0.1 
Z<-x
tol_error<-0.1
max_iter<-100

sc_bicluster <- function(lambda,gamma1,gamma2,del,Z,X,u,v,tol_error,max_iter){
  
  #lambda/gamma2 = tunning parameter of SCAD penalty
  
  n<- nrow(X)
  p<- ncol(X)
  x_old <- X
  u_old <-u
  v_old<-v
  
  iter_error <-matrix(ncol=5,nrow=max_iter)%>%
    'colnames<-'(value=c("X","A","B","U","V"))
  
  
  # initialize A_ij (X_i-X_j (ith row of X))
  xt<-t(x_old)
  vecxt<-as.matrix(c(xt))
  A<-matrix(0,nrow=p,ncol=n*(n-1)*(1/2))
  E_ij<-list()
  m1<-combinations(n,2)
  for (idx in 1:nrow(m1)){
    i=m1[idx,1]
    j=m1[idx,2]
    Ei=matrix(0,nrow=p,ncol=n*p)
    Ej=matrix(0,nrow=p,ncol=n*p)
    for(k in 1:p){
      Ei[k,(i-1)*p+k]<-1
      Ej[k,(j-1)*p+k]<-1
    }
    E_ij[[idx]]=Ei-Ej
    A[,idx]<-E_ij[[idx]]%*%vecxt
    
    
  }
  A_old<-A
  
  
  # initialize B_ij (X_i-X_j (ith column of X))
  vecx=c(x_old)
  B<- matrix(0,nrow=n,ncol=p*(p-1)*(1/2))
  m2<-combinations(p,2)
  E_kl<-list()
  for (idx in 1:nrow(m2)){
    i=m2[idx,1]
    j=m2[idx,2]
    Ei=matrix(0,nrow=n,ncol=n*p)
    Ej=matrix(0,nrow=n,ncol=n*p)
    for(k in 1:n){
      Ei[k,(i-1)*n+k]<-1
      Ej[k,(j-1)*n+k]<-1
      }
    E_kl[[idx]]=Ei-Ej
    B[,idx]<-E_kl[[idx]]%*%vecx

    
  }
  B_old<-B
  
  
  
  #index matrix of X_omega
  idxom<-matrix(0,nrow=n*p,ncol=n*p)
  #m3 <- expand.grid(1:n,1:p)
  vecz<-c(Z)
  for (idx in 1:length(vec(z))){
      if(vecz[idx]==0)
        diag(idxom)[idx]=1
      else
        diag(idxom)[idx]=0
      
    }
  
  for(iter in 1:max_iter){
  
    #update X
    
    
    summat1=matrix(0,nrow=n*p,ncol=n*p)
    summat3=matrix(0,nrow=n*p,ncol=1)
    summat4=matrix(0,nrow=n*p,ncol=1)
    for (idx in 1:nrow(m1)){
      summat1=summat1+t(E_ij[[idx]])%*%(E_ij[[idx]])
      summat3=summat3+t(E_ij[[idx]])%*%u_old[,idx]
      summat4=summat4+t(E_ij[[idx]])%*%A_old[,idx]
    }
    
    
    
    summat2=matrix(0,nrow=n*p,ncol=n*p)
    summat5=matrix(0,nrow=n*p,ncol=1)
    summat6=matrix(0,nrow=n*p,ncol=1)
    
    for (idx in 1:nrow(m2)){
  
      summat2=summat2+t(E_kl[[idx]])%*%E_kl[[idx]]
      summat5=summat5+t(E_kl[[idx]])%*%v_old[,idx]
      summat6=summat5+t(E_kl[[idx]])%*%B_old[,idx]
    }
    
    
    
    vecx_new=solve(diag(n*p)-gamma1*t(idxom)%*%idxom+del*summat1+del*summat2)%*%(c(Z)-summat3+del*summat4-summat5+del*summat6)
    
    #update A
    
    A_new<-data.frame(matrix(0,nrow=p,ncol=n*(n-1)*(1/2)))
    norm_vec <- function(x) sqrt(sum(x^2))
    for (idx in 1:nrow(m1)){
      gamma_ij=norm_vec(E_ij[[idx]]%*%vecx_new-(1/del)*u_old[,idx])
      a_ij=del*E_ij[[idx]]%*%vecx_new-u_old[,idx]
      
      if(norm_vec(A_old[,idx])<=lambda){
        if(-lambda<(del*E_ij[[idx]]-u_old[,idx])&&(del*E_ij[[idx]]-u_old[,idx])<lambda){
          A_new[,idx]=0
        }
        else{
          A_new[,idx]=(1-(lambda/(2*del*gamma_ij)))*(E_ij[[idx]]%*%vecx_new-(1/del)*u_old[,idx])
          
        }
        
      }
      
      else if(lambda<norm_vec(E_ij[[idx]])&&norm_vec(E_ij[[idx]]<gamma2*lambda)){
    
          A_new[,idx]=((1/del)-(gamma2*lambda/(2*del*(gamma2-1)*norm_vec(a_ij))))*a_ij
         
      }
      else
        A_new[,idx]=E_ij[[idx]]%*%vecx_new-(1/del)*u_old[,idx]
      
        
      
    }
    
    
    #update B
    
    B_new<-data.frame(matrix(0,nrow=n,ncol=p*(p-1)*(1/2)))
    for (idx in 1:nrow(m2)){
      gamma_kl=norm_vec(E_kl[[idx]]%*%vecx_new-(1/del)*v_old[,idx])
      a_kl=del*E_kl[[idx]]%*%vecx_new-v_old[,idx]
      if(norm_vec(B_old[idx])<=lambda){
        if(-lambda<(del*E_kl[[idx]]-v_old[,idx])&&(del*E_kl[[idx]]-v_old[,idx])<lambda){
          B_new[,idx]=0
        }
        else{
          
          
          B_new[,idx]=(1-(lambda/(2*del*gamma_kl)))*(E_kl[[idx]]%*%vecx_new-(1/del)*v_old[,idx])
          
        }
        
      }
      
      else if(lambda<norm_vec(E_kl[[idx]])&&norm_vec(E_kl[[idx]]<gamma2*lambda)){
        
          B_new[,idx]=((1/del)-(gamma2*lambda/(2*del*(gamma2-1)*norm_vec(a_kl))))*a_kl
      }
      else
        B_new[,idx]=E_kl[[idx]]%*%vecx_new-(1/del)*v_old[,idx]
      
      
      
    }
    
    #update u,v
    
    u_new=matrix(0,nrow=p,ncol=n*(n-1)*(1/2))
    v_new=matrix(0,nrow=n,ncol=p*(p-1)*(1/2))
    
    for (idx in 1:nrow(m1)){
      u_new[,idx]=u_old[,idx]+del*(A_new[,idx]-E_ij[[idx]]%*%vecx_new)
    }
    for(idx in 1:nrow(m2)){
      v_new[,idx]=v_old[,idx]+del*(B_new[,idx]-E_kl[[idx]]%*%vecx_new)
    }
    
    
    #iteration error
    x_new<-matrix(vecx_new,nrow=n,ncol=p)
    iter_error[iter,"X"] <-Matrix:: norm(x_old-x_new,type="F")
    iter_error[iter,"A"] <-Matrix::norm(as.matrix(A_old-A_new),type="F") 
    iter_error[iter,"B"] <-Matrix::norm(as.matrix(B_old-B_new),type="F")
    iter_error[iter,"U"] <-Matrix:: norm(as.matrix(u_old-u_new),type="F")
    iter_error[iter,"V"] <-Matrix:: norm(as.matrix(v_old-v_new),type="F")
    
    if(sum(iter_error[iter,])<tol_error) break
    
    
    
    x_old<-x_new
    A_old<-A_new
    B_old<-B_new
    u_old<-u_new
    v_old<-v_new
  }
  
  return(list(X=x_new,A=A_new,B=B_new,U=u_new,V=v_old,iter_error=iter_error))

}
