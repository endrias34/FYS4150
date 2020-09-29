#ifndef INITIALIZE_H
#define INITIALIZE_H

void fill_beam_array(arma::mat& A, arma::mat& C, arma::vec& B, int n);
void fill_one_electron_array(arma::mat& A, int n, int rhomax);
void fill_two_electrons_array(arma::mat& A, int n, double rhomax, double omega);

#endif // INITIALIZE_H
