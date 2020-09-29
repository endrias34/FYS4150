# include <iostream>
# include <cmath>
# include <armadillo>
# include "jacobi.h"

/* Testing max a(i,j)*/
void test_max_offdiag()
{
    int k, l, N = 5;
    double max_val = -17.4;              // Choose a max value
    double max_cal;
    arma::mat A = arma::zeros(N, N);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A(i,j) = i*j*j - i;    // Some arbitrary values
            }

            else {
                if (i*j == 12){
                    A(i,j) = max_val;  // Fill max value here
                }

                else {
                    A(i,j) = i*j - j;   // Some arbitrary values
                }
            }

            A(j,i) = A(i, j);  // Make sure matrix is symmetric
        }
    }
    max_cal = max_offdiag(A, k, l, N);
    if (max_cal == fabs(max_val)){
        std::cout << "Max non-diagonal element found ====> TEST PASSED!" << std::endl;
    }

    else {
        std::cout << "Max non-diagonal element NOT found ==== > TEST FAILED!" << std::endl;
    }
}

/* Testing eigenvalues */
void test_jacobi_eigen()
{
    int N        = 3;
    int niter    = 0;
    double eig_0 = -0.113538;       
    double eig_1 = 1.59642;  
    double eig_2 = 5.51711;   
    double eps   = 1.0E-5;       

    arma::mat A = arma::zeros(N, N);    
    arma::mat V = arma::eye(N, N);    

    A(0,0) = 2.0;  A(1,0) = 1.0;  A(2,0) = -1.0;  
    A(0,1) = 1.0;  A(1,1) = 4.0;  A(2,1) = -2.0;
    A(0,2) = -1.0; A(1,2) = -2.0; A(2,2) = 1.0;
 
    jacobi_eigen(A, V, N, niter);    
    arma::vec eigval_jac = arma::sort(A.diag());  

    if (eigval_jac(0) - eig_0 > eps){
        std::cout << "Eigenvalue 1 is wrong ====> TEST FAILED!" << std::endl;
    }

    else if (eigval_jac(1) - eig_1 > eps){
        std::cout << "Eigenvalue 2 is wrong ====> TEST FAILED!" << std::endl;
 
    }

    else if (eigval_jac(2)- eig_2 > eps){
        std::cout << "Eigenvalue 3 is wrong ====> TEST FAILED!" << std::endl;
    }

    else {
        std::cout << "Correct eigenvalues found =========> TEST PASSED!" << std::endl;
    }
}


/* Testing eigenvector orthogonality */
void test_orthogonality(){
    int N     = 5;
    int niter = 0;
    arma::mat A = arma::zeros(N, N);    // Matrix to solve for
    arma::mat V = arma::eye(N, N);      // Matrix to hold eigenvectors

    /* Fill A with some values. */
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A(i,j) = i*j*j - i;    // Some arbitrary values
            }

            else {
                A(i,j) = i*j - j;   // Some arbitrary values
                A(j,i) = A(i, j);  // Make sure matrix is symmetric
            }
        }
    }

    jacobi_eigen(A, V, N, niter);    // Solve eigenvalue problem

    bool test_failed = false;
    double eps = 1.0E-10;       // Tolerance for "close-to-zero"
    double inner_prod = 0.0;


    for (int c_1 = 0; c_1 < N; c_1++){      // Loop through cols in V.
        for (int c_2 = 0; c_2 < N; c_2++){  // Loop through cols in V per col in V
            inner_prod = 0.0;

            for (int k = 0; k < N; k++){
                inner_prod += V(k,c_1)*V(k,c_2);    // Compute inne product
            }

            if (inner_prod > eps && c_1 != c_2){
                test_failed = true;
            }

            if (inner_prod - 1 > eps && c_1 == c_2){
                test_failed = true;
            }
        }
    }

    if (test_failed){
        std::cout << "Orthogonality is NOT preserved ========> TEST FAILED" << std::endl;
    }

    else {
        std::cout << "Orthogonality is preserved ========> TEST PASSED!" << std::endl;
    }
}
