
#ifndef SOLVER_H
#define SOLVER_H

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <cmath>
# include <ctime>
# include "planet.h"
# include "All_const.h"

class solver
{
    public:

        solver();
        solver(planet* n_planets, int n, bool state_sun = false);
        void simulate(double h, double t_max, double t_write, std::string method, std::string casenr);
        void write_to_file(int file_index, double t, double x, double y, double z);
        void center_of_mass() const;
        ~solver();

    private:

        void euler(double h, double t_max, int frame_write);
        void eulerTime(double h, double t_max);
        void eulerEnergy(double h, double t_max, int frame_write);
        void eulerAngmom(double h, double t_max, int frame_write);
        void verlet(double h, double t_max, int frame_write);
        void verletTime(double h, double t_max);
        void verletEnergy(double h, double t_max, int frame_write);
        void verletAngmom(double h, double t_max, int frame_write);
        void verlet_GR(double h, double t_max, int frame_write);
        double total_acc(planet subject, planet* objects, int dim) const;
        double total_acc_GR(planet subject, int dim);
        void energy(double time,int frame_write, int frame) const;
        void angular_momentum(double time, int frame_write, int frame) const;
        

        // properties
        int num_planets;
        bool stationary_sun;
        planet* n_planets;
        std::ofstream* ofiles;
        

};

#endif // SOLVER_H
