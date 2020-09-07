# include <iostream>
# include <cmath>
# include <string>
# include <iomanip>
# include <fstream>
# include <cstdlib>
# include <ctime>
# include <armadillo>

/* ------------------------------------------------------------------------------------------------------------ */
/* Declare functions that are going to be used later */
void initialize(double *&, double *&, double *&, double *&, double *&, double *&, double *&, double *&, double *&, int);
void sol_general(double*, double*, double*, double*, int, double*);
void sol_special(double*, double*, int, double*);
void sol_lu(int);
double sol_analytic(double);
double source_term(double);
void time_algorithm(std::string, int, int);
void write_to_file(const char*, double*, double*, int);
/* ------------------------------------------------------------------------------------------------------------ */

int main(int argc, char* argv[])
{
    for (int i = 1; i < argc; i++)
    {
        int n = atoi(argv[i]);  // Number of grid points

        std::string fileoutg = "Result_"; // output filename definition
        fileoutg.append(argv[i]);            
        fileoutg.append("_general.txt");

        std::string fileouts = "Result_"; 
        fileouts.append(argv[i]);            
        fileouts.append("_special.txt");

        double *a, *bg, *bs, *c, *x, *gg, *gs, *num_solg, *num_sols;    // Pointer variable declaration.

        initialize(a, bg, bs, c, x, gg, gs, num_solg, num_sols, n);  // memory allocation and initialization

        /* Solving traidaigonal linear system of eqations  and write to file */
        sol_general(a, bg, c, gg, n, num_solg);       // general traidiagonal matrix algorithm
        sol_special(bs, gs, n, num_sols);             // specialized (Topliz) traidiagonal matrix algorithm
        write_to_file(fileoutg.c_str(), x, num_solg, n);    // write results to file
        write_to_file(fileouts.c_str(), x, num_sols, n);    

        sol_lu(n);   // LU-decomposition based algorithm

        delete[] a; delete[] bg; delete[] bs; delete[] c; delete[] x; delete[] gg; delete[] gs; 
        delete[] num_solg; delete[] num_sols;   // Free up memory
   
    }
    
    time_algorithm("general", 3, 1000);   // Run timing for general algorithm
    time_algorithm("special", 3, 1000);   // Run timing for special algorithm
    
    return 0;   
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Function that allocates memory and initialize diagonals a, b, c, x, g and the solution */
void initialize(double *& a, double *& bg, double *& bs, double *& c, double *& x, double *& gg, double *& gs, 
double *& solg, double *& sols, int n)
{
    a = new double[n-1];              
    bg = new double[n];                 
    bs = new double[n];                 
    c = new double[n-1];               
    x = new double[n];                
    gg = new double[n];                 
    gs = new double[n];                 
    solg = new double[n];              
    sols = new double[n];             

    double h = 1.0/((double)(n + 1));  
    double h_squared = h*h;            

    for (int i = 0; i < n-1; i++) a[i] = -1.0;     
    for (int i = 0; i < n-1; i++) c[i] = -1.0;     
    for (int i = 0; i < n; i++) bg[i] = 2.0;         
    for (int i = 0; i < n; i++) bs[i] = 2.0;        
    for (int i = 0; i < n; i++) x[i] = (i+1)*h;     
    for (int i = 0; i < n; i++) gg[i] = h_squared*source_term(x[i]);   
    for (int i = 0; i < n; i++) gs[i] = h_squared*source_term(x[i]);   
}
/* ------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------ */
/* General traidiagonal matrix algorithm */
void sol_general(double* a, double* bg, double* c, double* yg, int n, double* solutiong)
{
    
    /* Forward substitution. */
    for (int i = 1; i < n; i++){
        bg[i] = bg[i] - a[i-1]*c[i-1]/bg[i-1];     
        yg[i] = yg[i] - (a[i-1]/bg[i-1])*yg[i-1];   
    }

    /* Backward substitution. */
    solutiong[n-1] = yg[n-1]/bg[n-1];  

    for (int i = n-2; i >= 0; i--){
        solutiong[i] = (yg[i] - c[i]*solutiong[i+1])/bg[i];  
    }
}
/* ------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------ */
/*  Special traidiagonal matrix algorithm */
void sol_special(double* bs, double* ys, int n, double* solutions)
{
    /* Forward substitution. */
    for (int i = 1; i < n; i++){
        bs[i] = (i + 2)/((double)(i + 1));       
        ys[i] = ys[i] + (ys[i-1]/bs[i-1]);      
    }

    /* Backward substitution. */
    solutions[n-1] = ys[n-1]/bs[n-1];          

    for (int i = n-2; i >= 0; i--){
        solutions[i] = (ys[i] + solutions[i+1])/bs[i];  
    }
}
/* ------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------ */
/* LU-decomposition based algorithm using Armadillo */
void sol_lu(int n)
{
    std::ofstream outfilelu;                  
    std::string time_filenamelu = "timing_";  
    time_filenamelu.append("LU");        
    time_filenamelu.append(".txt");
    outfilelu.open(time_filenamelu);
    outfilelu << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    outfilelu << std::setw(8) << "N" << std::setw(20) << "time" << std::endl;      

    int upper_limit = 4;     

    for (int k = 1; k < upper_limit; k++)
    {
        int N = (int)pow(10.0, (double)(k)); 
        double h = 1.0/((double)(N+1));  
        double h_squared = h*h;         
        
        std::string s = std::to_string(N);
        std::string fileoutlu = "Result_"; 
        fileoutlu.append(s);            
        fileoutlu.append("_LU.txt");        

        arma::mat A(N, N);      
        arma::vec x(N);         
        arma::vec g(N);         

        /* Initialize x ang g */
        for (int i = 0; i < N; i++) x[i] = (i+1)*h;
        for (int i = 0; i < N; i++) g[i] = h_squared*source_term(x[i]);

        /* Initialize matrix A */
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                if (i == j){
                    A(i, j) = 2.0;              
                }

                else if (fabs(i - j) == 1){
                    A(i, j) = -1.0;             
                }

                else {
                    A(i, j) = 0.0;      
                }
            }
        }

        arma::vec v(N);                

        double num_runs_per_N = 1000;
        double time_total = 0.0;
        clock_t start_time, end_time;

        /* Run Armadillo several times to get representative average runtime. */
        for (int i = 0; i < num_runs_per_N; i++)
        {
            start_time = clock();           
            v = arma::solve(A, g);          
            end_time = clock();            
            double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
            time_total += time_used;
        }

        double time_averagelu = time_total/num_runs_per_N;

        outfilelu << std::setw(10) << N;
        outfilelu << std::setw(20) << std::setprecision(8) << time_averagelu << std::endl;

        double* solutionlu = new double[N];
        for (int i = 0; i < N ; i++)
        {
            solutionlu[i] = v[i];              
        }
        if(n == N){
        write_to_file(fileoutlu.c_str(), x.memptr(), solutionlu, N);    // Write solution to file
        }
        delete[] solutionlu;
    }
}
/* ------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------ */
/* Function returning the value of the source term for a given x. */
double source_term(double x)
{
    return 100.0*exp(-10.0*x);
}

/* Function returning the ecaxt solution for a given. */
double sol_analytic(double x)
{
    return 1 - (1 - exp(-10.0))*x - exp(-10.0*x);
}
/* ------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------ */
/* Function that runs the requested algorithm for a series
 * of N-values and then writes the runtime to file. */
void time_algorithm(std::string algorithm, int num_Ns, int num_runs_per_N)
{
    std::ofstream outfile;                  
    std::string time_filename = "timing_";  
    time_filename.append(algorithm);        
    time_filename.append(".txt");
    outfile.open(time_filename);
    outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    outfile << std::setw(8) << "N" << std::setw(20) << "time" << std::endl;  

    for (int i = 0; i < num_Ns; i++){
        int N = (int)pow(10.0, (double)(i+1));       
        double time_total = 0.0;                   

        for (int j = 0; j < num_runs_per_N; j++){
            double *a, *bg, *bs, *c, *x, *gg, *gs, *solutiong, *solutions;      
            initialize(a, bg, bs, c, x, gg, gs, solutiong, solutions, N);     

            clock_t start_time = clock();          

            if (!algorithm.compare("general")){
                sol_general(a, bg, c, gg, N, solutiong);  
            }

            else if (!algorithm.compare("special")){
                sol_special(bs, gs, N, solutions);       
            }

            clock_t end_time = clock();            
            double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
            time_total += time_used;   

            delete[] a; delete[] bg; delete[] bs; delete[] c; delete[] x; delete[] gg; delete[] gs;
            delete[] solutiong; delete[] solutions; 
        }

        double time_average = time_total/num_runs_per_N;    
        
        /* Write time data to file. */
        outfile << std::setw(10) << N;
        outfile << std::setw(20) << std::setprecision(8) << time_average << std::endl;
    }

    outfile.close();
}

/* ------------------------------------------------------------------------------------------------------------ */
void write_to_file(const char* filename, double* x, double* numerical, int N){
    std::ofstream ofile;    
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "              x:             Numerical_sol:      Analytical_sol:    Relative error" << std::endl;

    for (int i = 0; i < N; i++) {
        double x_val = x[i];
        double exact_val = sol_analytic(x[i]);
        double rel_error = fabs((exact_val - numerical[i])/exact_val);
        ofile << std::setw(20) << std::setprecision(8) << x_val;
        ofile << std::setw(20) << std::setprecision(8) << numerical[i];
        ofile << std::setw(20) << std::setprecision(8) << exact_val;
        ofile << std::setw(20) << std::setprecision(8) << rel_error << std::endl;
    }

    ofile.close();
}