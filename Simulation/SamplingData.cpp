// [[Rcpp::depends(Rcpp, RcppArmadillo, BH)]]

#include <RcppArmadillo.h>
#include <boost/math/special_functions/bessel.hpp>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
//Eigenvalues for a dense square matrix
arma::vec Eigenvalues(arma::mat X) {
  
  arma::cx_vec eigval;
  arma::eig_gen(eigval, X);
  
  return arma::conv_to<arma::vec>::from(eigval);
}

// [[Rcpp::export]]
//Eigenvectors for dense general square nonsymmetric matrix
arma::mat Eigenvectors(arma::mat X) {
  
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, X);
  
  arma::mat eigvc = arma::conv_to<arma::mat>::from(eigvec);
  return  -1 * eigvc;
}

// [[Rcpp::export]]
//Eigenvectors for a dense symmetric matrix
arma::mat EigenvectorsSym(arma::mat X) {
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X, "dc");
  
  arma::mat eigvc = arma::conv_to<arma::mat>::from(eigvec);
  return  -1 * eigvc;
}

// [[Rcpp::export]]
//Check if a matrix is identity
//Note that this function only works for square matrices
int IsIdentity(arma::mat X){
  //count rows and columns in X
  int ncols = X.n_cols;
  int nrows = X.n_rows;
  
  //create a testvector 
  
  arma::vec testI(ncols*nrows);
  
  //set a count
  int count = 0;
  
  //start checking
  for(int iii = 0; iii < ncols; ++iii){
    for(int jjj = 0; jjj < nrows; ++jjj){
      int testsim = X(jjj,iii);
      int testsim2 = X(iii,jjj);
      if(((testsim == 1) & (iii == jjj)) | ((testsim2 == 0) & (iii != jjj))){
        testI(count) = 0;
      }else{
        testI(count) = 1;
      }
      count = count + 1;
    }
  }
  //return sum of the test vector 
  //if sum of this vector is zero, the matrix is an identity matrix
  return sum(testI);
}

// [[Rcpp::export]]
//Sample from a multivariate normal distribution
//Note that sigma has to be a square matrix
arma::mat mvrnormArmaEigen(int n, arma::vec mu, arma::mat sigma){
  int ncols = sigma.n_cols;
  arma::mat Y(n, ncols);
  for(int iii=0; iii < ncols; ++iii){
    Y.col(iii) = as<arma::vec>(Rcpp::rnorm(n));
  }
  
  int testsym = 0;
  
  //test if matrix is idenity matrix
  if(IsIdentity(sigma) == 0){
    testsym = 1;
  }
  
  //test is matrix is symmetrical
  arma::mat sigmat = sigma.t();
  for(int iii=0; iii < ncols*ncols; ++iii){
    if(sigma(iii) != sigmat(iii)){
      testsym = 1;
      break;
    }
  }
  
  if(testsym == 1){
    return arma::repmat(mu, 1, n) + (EigenvectorsSym(sigma) * arma::diagmat(sqrt(max(Eigenvalues(sigma), arma::zeros(ncols)))) * Y.t());
  }else{
    return arma::repmat(mu, 1, n) + (Eigenvectors(sigma) * arma::diagmat(sqrt(max(Eigenvalues(sigma), arma::zeros(ncols)))) * Y.t());
  }
}


//[[Rcpp::export]]
//Sampling projected normal regression data with p predictors. Predictors are 
//centered after sampling outcome Y/Theta.
List RData(int N, int p, arma::vec Beta1, arma::vec Beta2, arma::vec mean, arma::vec sd){
  
  //Initialize matrices and vectors
  
  //Centered predictor matrices for both bivariate components.
  arma::mat X1(N, p);
  arma::mat X2(N, p);

  arma::mat Y(N,2);
  arma::vec Theta(N);
  arma::vec xc(N);
  arma::vec x(N);
  
  //Intercept
  arma::vec i = arma::ones<arma::vec>(N);

  X1.col(0) = i;
  X2.col(0) = i;
  
  //Sample predictors
  for(int kkk = 1; kkk < p; ++kkk){
    x = Rcpp::rnorm(N, mean(kkk-1), sd(kkk-1));
    X1.col(kkk) = x;
    X2.col(kkk) = x;
  }
  
  for(int iii = 0; iii < N; ++iii){
    
    arma::vec mu(2);
    arma::mat x1 = X1.row(iii) * Beta1;
    arma::mat x2 = X2.row(iii) * Beta2;
    mu(0) = x1(0,0);
    mu(1) = x2(0,0);
    Y.row(iii) = mvrnormArmaEigen(1, mu, arma::eye(2,2)).t();
    double y1 = Y(iii,0);
    double y2 = Y(iii,1);
    Theta.row(iii) = atan2(y2,y1);
  }
  
  return  Rcpp::List::create(Rcpp::Named("Theta") = Theta,
                             Rcpp::Named("X1") = X1,
                             Rcpp::Named("X2") = X2,
                             Rcpp::Named("Y") = Y);
}
