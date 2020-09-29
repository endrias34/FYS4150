
# include <iostream>
# include <armadillo>
# include <iomanip>
# include <fstream>
# include <string>
# include <ctime>
# include "unit_tests.h"
# include "jacobi.h"

int main(int argc, char* argv[])
{
    std::cout << "\nEXECUTING UNIT TESTS..." << std::endl
              << "==================================================" << std::endl;

    test_max_offdiag();
    test_jacobi_eigen();
    test_orthogonality();

    std::cout << "==================================================\n" << std::endl;

    return 0;
}

