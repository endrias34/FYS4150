# include <iostream>
#include <cmath>
# include <armadillo>


/* Function that initializes matrix for single electron case. */
void fill_beam_array(arma::mat& A, arma::mat& C, arma::vec& B, int n)
{  
  int i, j, Dim;
  double RMin, RMax, Step, DiagConst, NondiagConst; 
  double pi = acos(-1.0);
  RMin = 0.0; 
  RMax = 1.0; 
  Dim  = n;  
  
  // Integration step length
  Step         = RMax/ (Dim+1);
  DiagConst    = 2.0 / (Step*Step);
  NondiagConst = -1.0 / (Step*Step);
  
  // Analytical Eigen values
  for(i = 1; i <= n; i++) 
    {
      B(i-1)  = DiagConst + 2*NondiagConst*cos(i*pi*Step);
    
    }
 // Analytical Eigen vectors
  for(i = 1; i <= n; i++) {
      for(j = 1; j <= n; j++) {
          C(i-1,j-1) = sin(j*i*pi*Step);
        }
    }

  // Setting up tridiagonal matrix 
  A(0,0) = DiagConst;
  A(0,1) = NondiagConst;

  for(i = 1; i < Dim-1; i++) 
  {
    A(i,i-1)  = NondiagConst;
    A(i,i)    = DiagConst;
    A(i,i+1)  = NondiagConst;
  }

  A(Dim-1,Dim-2) = NondiagConst;
  A(Dim-1,Dim-1) = DiagConst;
}

/* Function that initializes matrix for single electron case. */
void fill_one_electron_array(arma::mat& A, int n, int rhomax)
{
    double rho_0 = 0.0;
    double rho_n = rhomax;
    double Step, DiagConst, NondiagConst;
    arma::vec rho(n);     

    // Integration step length
    Step         = rho_n/ (n+1);
    DiagConst    = 2.0 / (Step*Step);
    NondiagConst = -1.0 / (Step*Step);
 
    for (int i=0; i<n; i++)
    {
        rho(i) = rho_0 + (i+1)*Step;
    }

    // Setting up tridiagonal matrix 
    A(0,0) = DiagConst;
    A(0,1) = NondiagConst;

   for(int i = 1; i < n-1; i++) 
   {
     A(i,i-1)  = NondiagConst;
     A(i,i)    = DiagConst + (rho(i)*rho(i));
     A(i,i+1)  = NondiagConst;
   }

    A(n-1,n-2) = NondiagConst;
    A(n-1,n-1) = DiagConst + (rho(n-1)*rho(n-1));

}


/* Function that initializes matrix for two electrons case. */
void fill_two_electrons_array(arma::mat& A, int n, double rhomax, double omega)
{
    double rho_0 = 0.0;
    double rho_n = rhomax;
    double Step, DiagConst, NondiagConst;
    arma::vec rho(n);     

    // Integration step length
    Step         = (rho_n - rho_0)/ (n+1);
    DiagConst    = 2.0 / (Step*Step);
    NondiagConst = -1.0 / (Step*Step);
 
    for (int i=0; i<n; i++)
    {
        rho(i) = rho_0 + (i+1)*Step;
    }

    // Defining the frequency, 0.01, 0.5, 1, 5
    double omega_squared = omega*omega;

    // Setting up tridiagonal matrix 
    A(0,0) = DiagConst;
    A(0,1) = NondiagConst;

   for(int i = 1; i < n-1; i++) 
   {
     A(i,i-1)  = NondiagConst;
     A(i,i)    = DiagConst + (rho(i)*rho(i))*omega_squared + (1.0/rho(i));
     A(i,i+1)  = NondiagConst;
   }

    A(n-1,n-2) = NondiagConst;
    A(n-1,n-1) = DiagConst + (rho(n-1)*rho(n-1))*omega_squared + (1.0/rho(n-1));

}