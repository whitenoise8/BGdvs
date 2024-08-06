// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>  
#include <cmath>

using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat invtridiag(arma::mat A) {
  int n = A.n_cols;
  
  arma::vec alpha = zeros(n);
  arma::vec gamma = zeros(n);
  
  alpha(0) = A(0,0);
  
  double a21 = A(1,0);
  for (int i = 1; i < n; i++) {
    gamma(i) = a21/alpha(i-1);
    alpha(i) = A(i,i)-a21*gamma(i);
  }
  
  arma::mat C = zeros(n,n);
  C(n-1,n-1) = 1/alpha(n-1);  
  for (int j = n-2; j > -1; j--) {
    C(n-1,j) = -a21/alpha(j)*C(n-1,j+1);
    C(j,n-1) = C(n-1,j);
  }
  for (int i = n-2; i > 0; i--) for (int j = i-1; j > -1; j--) {
    C(i,i) = 1/alpha(i)+pow(a21/alpha(i),2)*C(i+1,i+1);
    C(i,j) = -a21/alpha(j)*C(i,j+1);
    C(j,i) = C(i,j);
  }
  C(0,0) = 1/alpha(0)+pow(a21/alpha(0),2)*C(1,1);
  
  return C;
}

// [[Rcpp::export]]
arma::mat invtridiag_memo(arma::vec dvec, double od, int n) {
  
  arma::vec alpha = zeros(n);
  arma::vec gamma = zeros(n);
  
  alpha(0) = dvec(0);
  
  double a21 = od;
  for (int i = 1; i < n; i++) {
    gamma(i) = a21/alpha(i-1);
    alpha(i) = dvec(i)-a21*gamma(i);
  }
  
  arma::mat C = zeros(n,n);
  C(n-1,n-1) = 1/alpha(n-1);  
  for (int j = n-2; j > -1; j--) {
    C(n-1,j) = -a21/alpha(j)*C(n-1,j+1);
    C(j,n-1) = C(n-1,j);
  }
  for (int i = n-2; i > 0; i--) for (int j = i-1; j > -1; j--) {
    C(i,i) = 1/alpha(i)+pow(a21/alpha(i),2)*C(i+1,i+1);
    C(i,j) = -a21/alpha(j)*C(i,j+1);
    C(j,i) = C(i,j);
  }
  C(0,0) = 1/alpha(0)+pow(a21/alpha(0),2)*C(1,1);
  
  return C;
}


// [[Rcpp::export]]
Rcpp::List getQuantitiesCpp(arma::cube A, arma::mat B, double nb, double bl) {
  
  arma::cube R = zeros(bl,bl,nb);
  R.slice(0) = inv_sympd(A.slice(0))*B;
  for (int j = 1; j < nb-1; j++) {
    R.slice(j) = inv(A.slice(j)-B*R.slice(j-1))*B;
  }
  
  arma::cube S = zeros(bl,bl,nb);
  S.slice(nb-2) = B*inv_sympd(A.slice(nb-1));
  for (int j = nb-3; j > -1; j--) {
    S.slice(j) = B*inv(A.slice(j+1)-S.slice(j+1)*B);
  }
  
  arma::cube D = zeros(bl,bl,nb);
  D.slice(0) = inv(A.slice(0)-B*S.slice(0).t());
  for (int j = 0; j < nb-2; j++) {
    D.slice(j+1) = inv(A.slice(j+1)-B*S.slice(j+1).t())*(eye(bl,bl)+B*D.slice(j)*S.slice(j));
  }
  D.slice(nb-1) = inv(A.slice(nb-1))*(eye(bl,bl)+B*D.slice(nb-2)*S.slice(nb-2));
  
  arma::mat oD = zeros(bl,nb);
  for (int i = 1; i < nb; i++) {
    arma::mat RD = R.slice(i-1)*D.slice(i);
    oD.col(i) = RD.diag();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("D") = D,
    Rcpp::Named("oD") = oD
  );
  
}

// [[Rcpp::export]]
Rcpp::List getXX(arma::mat mu_q_gamma, arma::mat X, arma::vec mu_q_s2inv) {
  double p = X.n_cols;
  double n = mu_q_s2inv.n_elem;
  
  Rcpp::List A(n);
  A[0] = zeros(p,p);
  for (int i = 1; i < n; i++) {
    arma::mat W = mu_q_gamma.row(i-1).t()*mu_q_gamma.row(i-1) + diagmat(mu_q_gamma.row(i-1)%(1-mu_q_gamma.row(i-1)));
    A[i] = X.row(i-1).t()*X.row(i-1)*mu_q_s2inv(i) % W;
  }
  
  return A;
}
