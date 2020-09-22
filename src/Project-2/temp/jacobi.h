#ifndef JACOBI_H
#define JACOBI_H

void max_offdiag(arma::mat &A,  int &k, int &l, int n);
double jacobi_rotation(arma::mat &A, arma::mat &R, int n, int k, int l);
void jacobi_eigen(arma::mat &A, arma::mat &R, int n, int &niter);

#endif // JACOBI_H