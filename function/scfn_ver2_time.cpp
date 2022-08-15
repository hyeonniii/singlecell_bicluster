#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <time.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List sc_bicluster2(double lambda, double gamma1, double del, arma::mat Z, arma::mat X, double tol_error, double max_iter){
  clock_t time_req=clock();
  
  double gamma2 = 3.7;
  int n= X.n_rows;
  int p= X.n_cols;
  arma::mat x_old=X;
  int m1=n*(n-1)/2;
  int m2=p*(p-1)/2;
  arma::mat u_old(p,m1, fill::ones);
  arma::mat v_old(n,m2, fill::ones);
  arma::vec vecx=vectorise(x_old);
  arma::mat n3(m1,2);
  arma::mat p3(m2,2);
  int k=0;
  for(int i=1; i<n;i++){
    for(int j=i+1;j<n+1;j++){
      
      
      n3(k,0)=i;
      
      n3(k,1)=j;
      k+=1;
    }
    
  }
  k=0;
  for(int i=1; i<p;i++){
    for(int j=i+1;j<p+1;j++){
      
      
      p3(k,0)=i;
      
      p3(k,1)=j;
      k+=1;
    }
    
  }
  
  
  List Eij(m1);
  for(int idx=0; idx<m1; idx++){
    int i=n3(idx,0);
    int j=n3(idx,1);
    IntegerVector k = seq(1,p);
    IntegerVector k1=(k-1)*n+i;
    IntegerVector k2=(k-1)*n+j;
    arma::mat E1(p,n*p);
    arma::mat E2(p,n*p);
    for(int i=0;i<p;i++){
      E1(k(i)-1,k1(i)-1)=1;
      E2(k(i)-1,k2(i)-1)=1;
    }
    Eij[idx]=E1-E2;
  }
  
  List Ekl(m2);
  for(int idx=0;idx<m2; idx++){
    int i=p3(idx,0);
    int j=p3(idx,1);
    IntegerVector k = seq(1,n);
    IntegerVector k1=(i-1)*n+k;
    IntegerVector k2=(j-1)*n+k;
    arma::mat E1(n,n*p);
    arma::mat E2(n,n*p);
    
    for(int i=0;i<n;i++){
      E1(k(i)-1,k1(i)-1)=1;
      E2(k(i)-1,k2(i)-1)=1;
    }
    Ekl[idx]=E1-E2;
  }
  
  clock_t time_req2=clock();
  clock_t time1=time_req2-time_req;
  Rcout<<"전처리(1)소요시간: "<<(float)time1/CLOCKS_PER_SEC <<" seconds\n";
  
  
  /* initialize A */
  
  
  
  
  arma::mat summat1(n*p,n*p);
  arma::mat summat2(n*p,n*p);
  arma::mat A_old(p,m1);
  arma::mat B_old(n,m2);
  for (int idx=0; idx<m1 ; idx++){
    arma::mat E=Eij[idx];
    arma::mat Et=trans(E);  
    summat1+=Et*E;
    A_old.col(idx)=E*vecx;
  }

  clock_t time_req3=clock();
  clock_t time2=time_req3-time_req2;
  Rcout<<"(2)initialize A: "<<(float)time2/CLOCKS_PER_SEC <<" seconds\n";
  /*initialize B*/
  
  
  for (int idx=0; idx<m2 ; idx++){
    arma::mat E=Ekl[idx];
    arma::mat Et=trans(E);
    summat2+=Et*E;
    B_old.col(idx)=E*vecx;
  }
  clock_t time_req4=clock();
  clock_t time3=time_req4-time_req3;
  Rcout<<"(3)initialize B: "<<(float)time3/CLOCKS_PER_SEC <<" seconds\n";

  
  /*index matrix of X_omega*/
  
  
  arma::vec vecz=vectorise(Z);
  arma::uvec col2=find(vecz);
  arma::mat idxom(n*p,n*p,fill::eye);
  for (uword i=0; i<col2.n_rows; i++){
    idxom(col2(i),col2(i))=0;
  }
  
  arma::mat iter_error(max_iter,5);
  
  arma::mat inv_mat=inv(id+2*gamma1*trans(idxom)*idxom+del*summat1+del*summat2);
  clock_t time_req5=clock();
  clock_t time4=time_req5-time_req4;
  Rcout<<"(4)initialize omega: "<<(float)time4/CLOCKS_PER_SEC <<" seconds\n";
  Rcout<< "finish initializing\n";
  
  
  /*iteration*/ 
  for(int iter=0;iter<max_iter;iter++){
    Rcout<<"iter : "<<iter+1<<"\n";
    Rcout<<"updating X\n";
    /*update X*/
    arma::mat summat3(n*p,1);
    arma::mat summat4(n*p,1);
    arma::mat summat5(n*p,1);
    arma::mat summat6(n*p,1);
    for (int i=0; i<m1;i++){
      arma::mat E=Eij[i];
      summat3+=trans(E)*u_old.col(i);
      summat4+=trans(E)*A_old.col(i);
    }
    for (int i=0; i<m2;i++){
      arma::mat E=Ekl[i];
      summat5+=trans(E)*v_old.col(i);
      summat6+=trans(E)*B_old.col(i);
    }
    arma::mat id(n*p,n*p,fill::eye);
    arma::vec vecx_new(n*p,1);
    vecx_new=inv_mat*(vecz+summat3+del*summat4+summat5+del*summat6);
    clock_t time_req6=clock();
    clock_t time5=time_req6-time_req5;
    Rcout<<"(5)update X: "<<(float)time5/CLOCKS_PER_SEC <<" seconds\n";
    
    /*update A*/
    Rcout<<"updating A\n";
    arma::mat A_new(p,m1);
    arma::vec a_ij(p);
    arma::vec A1(p);
    for(int i=0;i<m1;i++){
      a_ij=as<arma::mat>(wrap(Eij[i]))*vecx_new-(1/del)*u_old.col(i);
      if((norm(a_ij)>lambda/del)&(norm(a_ij)<=lambda*(del+1)/del)){
        A_new.col(i)=((del*norm(a_ij)-lambda)/(del*norm(a_ij)))*a_ij;
      }else if(norm(a_ij)<=lambda/del){
        A_new.col(i)=A1;
      }else if(norm(a_ij)<=gamma2*lambda){
        A_new.col(i)=((del*(gamma2-1)*norm(a_ij)-gamma2*lambda)/(norm(a_ij)*(del*(gamma2-1)-1)))*a_ij;
      }else{
        A_new.col(i)=a_ij;
      }
      
    }
    clock_t time_req7=clock();
    clock_t time6=time_req7-time_req6;
    Rcout<<"(6)update A: "<<(float)time6/CLOCKS_PER_SEC <<" seconds\n";
    
    /*update B*/
    Rcout<<"updating B\n";
    arma::mat B_new(n,m2);
    arma::vec b_ij(n);
    arma::vec B1(n);
    for(int i=0;i<m2;i++){
      b_ij=as<arma::mat>(wrap(Ekl[i]))*vecx_new-(1/del)*v_old.col(i);
      if((norm(b_ij)>lambda/del)&(norm(b_ij)<=lambda*(del+1)/del)){
        B_new.col(i)=((del*norm(b_ij)-lambda)/(del*norm(b_ij)))*b_ij;
      }else if(norm(b_ij)<=lambda/del){
        B_new.col(i)=B1;
      }else if(norm(b_ij)<=gamma2*lambda){
        B_new.col(i)=((del*(gamma2-1)*norm(b_ij)-gamma2*lambda)/(norm(b_ij)*(del*(gamma2-1)-1)))*b_ij;
      }else{
        B_new.col(i)=b_ij;
      }
      
    }
    clock_t time_req8=clock();
    clock_t time7=time_req8-time_req7;
    Rcout<<"(7)update B: "<<(float)time7/CLOCKS_PER_SEC <<" seconds\n";
    /*update U*/
    Rcout<<"updating U\n";
    arma::mat u_new(p,m1);
    arma::mat v_new(n,m2);
    for(int i=0;i<m1;i++){
      u_new.col(i)=u_old.col(i)+del*(A_new.col(i)-as<arma::mat>(wrap(Eij[i]))*vecx_new);
    }
    
    clock_t time_req9=clock();
    clock_t time8=time_req9-time_req8;
    Rcout<<"(8)update U: "<<(float)time8/CLOCKS_PER_SEC <<" seconds\n";
    
    Rcout<<"updating V\n";
    for(int i=0;i<m2;i++){
      v_new.col(i)=v_old.col(i)+del*(B_new.col(i)-as<arma::mat>(wrap(Ekl[i]))*vecx_new);
    }
    clock_t time_req10=clock();
    clock_t time9=time_req10-time_req9;
    Rcout<<"(9)update V: "<<(float)time9/CLOCKS_PER_SEC <<" seconds\n";
    
    
    arma::mat x_new=reshape(vecx_new,n,p);
    iter_error(iter,0)=norm(x_old-x_new,"fro");
    iter_error(iter,1)=norm(A_old-A_new,"fro");
    iter_error(iter,2)=norm(B_old-B_new,"fro");
    iter_error(iter,3)=norm(u_old-u_new,"fro");
    iter_error(iter,4)=norm(v_old-v_new,"fro");
    double error=sum(iter_error.row(iter));
    x_old=x_new;
    A_old=A_new;
    B_old=B_new;
    u_old=u_new;
    v_old=v_new;
    if(error<tol_error) break;
  }
  time_req=clock()-time_req;
  Rcout<<"소요시간: "<<(float)time_req/CLOCKS_PER_SEC <<" seconds";
  List result= List::create(Named("X")=x_old,Named("A")=A_old,Named("B")=B_old,Named("U")=u_old,Named("V")=v_old,Named("iter_error")=iter_error);
  return result;
}

