# include <iomanip>
# include <iostream>
# include <fstream>
# include <random>
# include <omp.h>
# include <algorithm>
# include <vector>




// g++ -fopenmp IsingSocialMain.cpp


// using namespace std;
std::random_device rdev;
std::mt19937_64 randng(rdev());
std::uniform_real_distribution<double> u_cont_d(0.0, 1.0);

//cd dropbox/fam/phd/fys4150/isingsocial/code

void vector(int*& Vptr, int n_agents);
void initialize(int n_agents, int order, double consentration, int* agents_system);
void MC(std::string filename, int mc_cycles, int n_agents, int* agents_system, double& M, int print_mc, double p_rules, int print_all_agents);

int main(int argc, char *argv[])
{
    int order;
    int n_timesteps;
    int systems_to_run;
    int n_agents;
    int mc_cycles;
    int print_mc = 0;
    int print_all_agents = 1;
    double p_rules = 1;
    double consentration = 0;
    std::string filename_start;

    n_agents = atoi(argv[1]);
    if (argc >= 3){order = atoi(argv[2]);}   // If input bigger than 0, you will print out results at each step you input
    if (argc >= 4){print_mc = atoi(argv[3]);}   // If input bigger than 0, you will print out results at each step you input
    if (argc >= 5){print_all_agents = atoi(argv[4]);}   // Input 1 to print out all spins, or input 0 to print out magnetization
    if (argc >= 6){consentration = atof(argv[5]);}  //  Input a concentration. 1 means all 1, and 0.9 means 90% 1's
    if (argc >= 7){p_rules = atof(argv[6]);}    // Probability of not following the neighbour rules, and an agent just randomly choose what to be


    mc_cycles = 100*n_agents;   // Found to give stable results


    int fileID = (int) (u_cont_d(randng)*(double)1000);
    if (print_mc == 0)
    {
        filename_start = "MResult";
        systems_to_run = (int) 1000 / omp_get_max_threads();
        omp_set_num_threads(omp_get_max_threads());
    }
    else
    {
        filename_start = "Result";
        systems_to_run = 1;
        omp_set_num_threads(1);
    }
    std::string filename = "Results/"+filename_start + std::to_string(n_agents) + "Agents" + std::to_string(mc_cycles) + "McCycles" + std::to_string(fileID) + ".txt";
    std::ifstream ifile;
    ifile.open(filename);



    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app); // append instead of overwrite


    // // Makes locks so only one process write to file
    omp_lock_t value_lock;
    omp_init_lock(&value_lock);


    #pragma omp parallel
    {
    for (int n = 0; n < systems_to_run; n++)
    {

        int *agents_system;
        double M = 0;
        vector(agents_system, n_agents);
        initialize(n_agents, order, consentration, agents_system);
        MC(filename, mc_cycles, n_agents, agents_system, M, print_mc, p_rules, print_all_agents);


        if (print_mc == 0)
        {
            for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
            {
                M += (int) agents_system[i];
            }
            omp_set_lock(&value_lock);
            outfile << M << "\n";
            omp_unset_lock(&value_lock);
        }
        delete [] agents_system;    // Release Memory
        std::cout << n << std::endl;

    }
    }
    outfile.close();

}







// function to initialize agent_vector
void vector(int*& Vptr, int n_agents)
{
    Vptr = (int *) malloc((n_agents+1)*sizeof(int));
}





// Initialize the agent_system
void initialize(int n_agents, int order, double consentration, int* agents_system)
{


    if (consentration == 0)
        {   // Initialize randomly from a uniform distribution
        for(int i = 0; i < (n_agents+1); i++) 
        {
            // order = 0 gives random ordering, and 1 gives all equal to -1
            agents_system[i] = (order == 1) ? 1 : (u_cont_d(randng) > 0.5) ? 1 : -1;
        }   // End of initialization
    }
    else
    {   // Initialize with concentration
        if (order == 0) // Randomized order
        {
            std::vector<int> v;
            for (int i = 0; i < n_agents; i++)
            {
                v.push_back(i);
            }
            std::shuffle(v.begin(), v.end(), randng);   // Make a shuffled index
            agents_system[0] = 0;
            for (int i = 0; i < n_agents; i++)
            {
                agents_system[v[i]+1] = (i+1 < consentration*n_agents) ? -1 : 1;
            }
        }   // End of randomizing order
        else
        {   // Initialize with a cluster
            agents_system[0] = 0;
            for (int i = 0; i < n_agents; i++)
            {
                agents_system[i+1] = (i+1 < consentration*n_agents) ? -1 : 1;
            }
        }   // End of cluster
    }   // End of concentration
}   // End Initialize








// Run MC-cycles
void MC(std::string filename, int mc_cycles, int n_agents, int* agents_system, double& M, int print_mc, double p_rules, int print_all_agents)
{
    int step_count = 0;

    std::ofstream outfile;
    if (print_mc != 0){outfile.open(filename, std::ios_base::app);} // append instead of overwrite

    if (p_rules == 1)
    {   // Start MC cylces
        for (int j = 0; j<mc_cycles; j++)
        {
            for (int i = 0; i<n_agents; i++)
            {

                step_count += 1;

                // Pick a random index
                int rand_i = (int) (u_cont_d(randng)*(double)(n_agents-1)) + 1 ;

                // Updating rules
                if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                {
                    agents_system[rand_i-1] = agents_system[rand_i];
                    agents_system[rand_i+2] = agents_system[rand_i];
                }
                if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                {
                    agents_system[rand_i-1] = agents_system[rand_i+1];
                    agents_system[rand_i+2] = agents_system[rand_i];
                }            
            }
        // Save magnetization or spins every step(= print_mc) cylce
        if (j % print_mc == 0)
        {

            if (print_all_agents == 0)
            {
                M = 0;
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    M += (double) agents_system[i];
                }
                outfile << M <<"\n";
            }
            else
            {
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    outfile << std::setw(3) << agents_system[i];
                }
                outfile << "\n";
            }
        }   // End of writing to file
        } // End of a cycle
    }   // End for loop of MC cycle



    else
    {   // Start MC cycles with a probability = p_rules to of not following rules
        for (int j = 0; j<mc_cycles; j++)
        {
            for (int i = 0; i<n_agents; i++)
            {

                step_count += 1;

                // Pick a random index
                int rand_i = (int) (u_cont_d(randng)*(double)(n_agents-1)) + 1 ;

                // Updating rules
                if (agents_system[rand_i]*agents_system[rand_i+1]==1)
                {
                    agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                    agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                }
                if (agents_system[rand_i]*agents_system[rand_i+1]==-1)
                {
                    agents_system[rand_i-1] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i+1] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                    agents_system[rand_i+2] = (u_cont_d(randng) > p_rules) ? agents_system[rand_i] : (u_cont_d(randng) > 0.5 ) ? 1 : -1;
                }            
            }
        // Save magnetization or spins every step(= print_mc) cylce
        if (j % print_mc == 0)
        {
            if (print_all_agents == 0)
            {
                M = 0;
                for(int i=1; i < n_agents+1; i++)   // Start with 1 since we initialize a system with index -1 = 0
                {
                    M += (double) agents_system[i];
                }
                outfile << M <<"\n";
            }
            else
            {
                for(int i=1; i < n_agents+1; i++)    // Start with 1 since we initialize a system with index -1 = 0
                {
                    outfile << std::setw(3) << agents_system[i];
                }
                outfile << "\n";
            }
        }   // End of writing to file

        }   // End of a cycle
    }   // End of MC cycles with probability of breaking the rules
}   // End MC cycles




