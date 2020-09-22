
# include <iostream>
# include <armadillo>
# include <iomanip>
# include <fstream>
# include <string>
# include <ctime>
# include "initialize.h"
# include "jacobi.h"

void write_results_to_file_plot(std::string fileout, arma::vec eig_vec_1, arma::vec eig_vec_2, arma::vec eig_vec_3, int n);

int main(int argc, char* argv[])
{
    
    std::string filename = argv[1];
    int n     = atoi(argv[2]);
    int niter = 0;

    std::string fileout = filename;

    arma::mat A  = arma::zeros(n,n);
    arma::mat R  = arma::eye(n,n);       
    arma::vec B  = arma::zeros(n);
    arma::vec C1 = arma::zeros(n);

    // Initialize matrix
    fill_beam_array(A,C1,B,n);

    // Time Armadillos Eigen solver.
    arma::vec eigval_ARMA;
    arma::mat eigvec_ARMA;
    clock_t start_time_ARMA = clock();
    arma::eig_sym(eigval_ARMA, eigvec_ARMA, A);
    clock_t end_time_ARMA = clock();
    double time_used_ARMA = (double)(end_time_ARMA - start_time_ARMA)/CLOCKS_PER_SEC;

    // Time jacobi implementation.
    clock_t start_time = clock();
    jacobi_eigen(A, R, n, niter);        
    clock_t end_time = clock();
    double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;

    std::cout << "Timing for " << n << "x" << n << " matrix:" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    std::cout << "Time used ARMADILLO: " << time_used_ARMA << std::endl;
    std::cout << "Time used Jacobi: " << time_used << std::endl;

    arma::vec eigval_JACOBI = arma::sort(A.diag());   // Sort eigenvalues in increasing order
    arma::uvec indices = arma::sort_index(A.diag());

    std::cout << "Jacobi rotaion nr.: " << niter << std::endl;
    
    //write_results_to_file_plot(fileout, B, eigval_JACOBI, eigval_ARMA, n);
    write_results_to_file_plot(fileout, C1, R.col(indices(0)), eigvec_ARMA.col(0), n);

    return 0;
}

void write_results_to_file_plot(std::string fileout, arma::vec eig_vec_1, arma::vec eig_vec_2, arma::vec eig_vec_3, int n)
{
    std::ofstream ofile;    // File object for output file
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "      Analytical:        Jacobi:           ARMADILLO:" << std::endl;
    for (int j = 0; j<n; j++){
        ofile << std::setw(30) << std::setprecision(20) << eig_vec_1(j);
        ofile << std::setw(30) << std::setprecision(20) << eig_vec_2(j);
        ofile << std::setw(30) << std::setprecision(20) << eig_vec_3(j) << std::endl;
     }

    ofile.close();
}
