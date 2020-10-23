
#ifndef PLANET_H
#define PLANET_H

# include <iostream>
# include <cmath>
# include "All_const.h"

class planet
{
    public:

        // Properties
        std::string name;
        double mass;
        double r[3], v[3], a[3];
        double beta;

        // Initializers
        planet();
        planet(std::string id, double m, double x, double y, double z, double v_x, double v_y, double v_z, double bb);
        
        // Functions
        double distance(planet planet_2) const;
        double acceleration(planet planet_2, int dim) const; 
        double kinetic_energy() const;
        double potential_energy(planet planet_2) const;
        double acc_mercury_GR(planet planet_2, int dim) const;

        ~planet();

};

#endif 
