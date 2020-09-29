# include <iostream>
# include <cmath>
# include <armadillo>

/* maximum value of the off diagonal elements. */
double max_offdiag(arma::mat A,  int &k, int &l, int n)
{
    double max_val = 0.0; 
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){      
            if (fabs(A(i,j)) > fabs(max_val)){
                max_val = fabs(A(i,j));
                k = i;                    
                l = j;       
            }
        }
    }

    return max_val;
}

/* Jacobi rotation algorithm*/
void jacobi_rotation(arma::mat &A, arma::mat &R, int n, int k, int l)
{
    double s, c;
    if ( A(k,l) != 0.0 ) 
    {
      double t, tau;
      tau = (A(l,l) - A(k,k))/(2*A(k,l));
    
      if ( tau >= 0 ) 
      {
        t = 1.0/(tau + sqrt(1.0 + tau*tau));
      } 
    
      else 
      {
        t = -1.0/(-tau +sqrt(1.0 + tau*tau));
      }
    
      c = 1/sqrt(1+t*t);
      s = c*t;
    } 

    else
    {
      c = 1.0;
      s = 0.0;
    }

    double a_kk = A(k,k);    
    double a_ll = A(l,l);    
    double a_ik, a_il, r_ik, r_il;

    A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
    A(k,l) = 0.0;    
    A(l,k) = 0.0;    

    for (int i = 0; i < n; i++)
    {
        if ((i != k) && (i != l))
        {
            a_ik   = A(i,k);       
            a_il   = A(i,l);       
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);      
            A(i,l) = a_il*c+ a_ik*s;
            A(l,i) = A(i,l);      
        }

        r_ik   = R(i,k);
        r_il   = R(i,l);
        R(i,k) = r_ik*c - r_il*s;
        R(i,l) = r_il*c + r_ik*s;
    }

}

/* Function that operates as a master function which oversees and executes the
 * jacobi rotation algorithm by calling on the other functions in this file. */
void jacobi_eigen(arma::mat &A, arma::mat &R, int n, int &niter)
{
    int k, l;                                          
    double curr_max_nd = 1.0E+10;      
    double eps         = 10.0E-10;   
    int curr_iter      = 0;
    int max_iter       = n*n*n;

    while ((curr_max_nd > eps) && (curr_iter < max_iter))
    {
        curr_max_nd = max_offdiag(A, k ,l, n);  
        jacobi_rotation(A, R, n, k, l);   
        niter = curr_iter;
        curr_iter++;
    }
}
