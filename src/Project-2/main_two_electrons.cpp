
# include <iostream>
# include <armadillo>
# include <iomanip>
# include <fstream>
# include <string>
# include <ctime>
# include "initialize.h"
# include "jacobi.h"

void write_to_file(std::string fileout, arma::vec out_2, arma::vec out_3, int n);
void write_to_file2(std::string fileout, arma::vec out_2, arma::vec out_3,
arma::vec out_5, arma::vec out_6, arma::vec out_8, arma::vec out_9, int n);

int main(int argc, char* argv[])
{
    std::string filename_evec = "/Users/endriasa/Projects/FYS4150/src/Project2/DataFiles/Result_eigenvec_2el_"; // output filename definition
    filename_evec.append(argv[1]);  
    filename_evec.append("_");
    filename_evec.append(argv[2]);    
    filename_evec.append("_");   
    filename_evec.append(argv[3]);        
    filename_evec.append(".txt");

    std::string filename_eval = "/Users/endriasa/Projects/FYS4150/src/Project2/DataFiles/Result_eigenval_2el_"; // output filename definition
    filename_eval.append(argv[1]); 
    filename_eval.append("_");
    filename_eval.append(argv[2]);    
    filename_eval.append("_");    
    filename_eval.append(argv[3]);      
    filename_eval.append(".txt");

    // Define nr. of grid points and nr. of Jacobi rotations
    int n         = atoi(argv[1]);
    double rhomax = atof(argv[2]);
    double omega  = atof(argv[3]);
    int niter     = 0;

    // Declare matrices and vectors
    arma::mat A = arma::zeros(n,n);
    arma::mat R = arma::eye(n,n);       
    arma::vec B = arma::zeros(n);
    arma::mat C = arma::zeros(n,n);

    // Initialize the system matrix
    fill_two_electrons_array(A, n, rhomax, omega);

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

    // Print the time used by Armadillo and Jacobi rotation
    std::cout << "N = "<<n <<" nr. of rotaions = " << niter<<" ARMADILLO Time = "<< time_used_ARMA << " Jacobi Time = " << time_used << " rho max = " << rhomax << " omega = " << omega << std::endl;

    // Sort eigenvalues in increasing order
    arma::vec eigval_JACOBI = arma::sort(A.diag());   
    arma::uvec indices      = arma::sort_index(A.diag());

    // Write to file the results
    write_to_file(filename_eval, eigval_JACOBI, eigval_ARMA, n);
    write_to_file2(filename_evec, R.col(indices(0)), eigvec_ARMA.col(0), R.col(indices(1)), eigvec_ARMA.col(1), 
    R.col(indices(2)), eigvec_ARMA.col(2), n);
 
    return 0;
}

void write_to_file(std::string fileout, arma::vec out_2, arma::vec out_3, int n)
{
    std::ofstream ofile;    
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "      Jacobi:               ARMADILLO:" << std::endl;
    for (int j = 0; j<n; j++){
        ofile << std::setw(30) << std::setprecision(20) << out_2(j);
        ofile << std::setw(30) << std::setprecision(20) << out_3(j) << std::endl;
     }

    ofile.close();
}
void write_to_file2(std::string fileout, arma::vec out_2, arma::vec out_3, arma::vec out_5, arma::vec out_6, 
arma::vec out_8, arma::vec out_9, int n)
{
    std::ofstream ofile;    
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "      Jacobi-1:               ARMADILLO-1:        Jacobi-2:               ARMADILLO-2:            Jacobi-3:               ARMADILLO-3:" << std::endl;
    for (int j = 0; j<n; j++){

        ofile << std::setw(30) << std::setprecision(20) << out_2(j);
        ofile << std::setw(30) << std::setprecision(20) << out_3(j); 
        ofile << std::setw(30) << std::setprecision(20) << out_5(j);
        ofile << std::setw(30) << std::setprecision(20) << out_6(j);
        ofile << std::setw(30) << std::setprecision(20) << out_8(j);
        ofile << std::setw(30) << std::setprecision(20) << out_9(j)<< std::endl;
     }

    ofile.close();
}  