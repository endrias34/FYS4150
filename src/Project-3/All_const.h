
#ifndef ALL_CONST_H
#define ALL_CONST_H

/* Namespace for all the physical and mathematical constants */
namespace consts 
{
    const double mass_sun     = 2.0E+30;
    const double pi           = acos(-1.0);
    const double four_pi_sq   = 4*acos(-1.0)*acos(-1.0);
    const int num_dims        = 3;                       
    const double G            = 6.67408E-11*2.97449433E-19;  
    const double c            = 173.0*365.0;   
    const double csq          = (173.0*365.0)*(173.0*365.0);          
    const double mercury_year = 90.0/365.0;     
    const double vel_esc      = 2*acos(-1.0)*sqrt(2);
}

#endif 
