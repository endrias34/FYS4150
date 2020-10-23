#include "planet.h"

planet::planet()
{
    this->name = "";
    this->mass = 0.0;
    this->r[0] = 0.0; this->r[1] = 0.0; this->r[2] = 0.0;
    this->v[0] = 0.0; this->v[1] = 0.0; this->v[2] = 0.0;
    this->a[0] = 0.0; this->a[1] = 0.0; this->a[2] = 0.0;
    this->beta = 2.0;
}

planet::planet(std::string id, double m, double x, double y,
               double z, double v_x, double v_y, double v_z, double bb)
{
    this->name = id;
    this->mass = m;
    this->beta = bb;
    this->r[0] = x; this->r[1] = y; this->r[2] = z;
    this->v[0] = v_x; this->v[1] = v_y; this->v[2] = v_z;
    this->a[0] = 0.0; this->a[1] = 0.0; this->a[2] = 0.0;

}

double planet::distance(planet planet_2) const
{
    double dist = pow(this->r[0] - planet_2.r[0], 2.0) 
                + pow(this->r[1] - planet_2.r[1], 2.0) 
                + pow(this->r[2] - planet_2.r[2], 2.0);

    return sqrt(dist);
}

double planet::acceleration(planet planet_2, int dim) const
{
    double r      = pow(this->distance(planet_2), this->beta+1.0);
    double mratio = planet_2.mass/consts::mass_sun;
    double delr   = this->r[dim] - planet_2.r[dim];
    double acc    = -(consts::four_pi_sq*mratio)*(delr/r);
    return acc;
}

double planet::acc_mercury_GR(planet planet_2, int dim) const
{
    double l_x = (this->r[1]*this->v[2] - this->r[2]*this->v[1]);
    double l_y = (this->r[2]*this->v[0] - this->r[0]*this->v[2]);
    double l_z = (this->r[0]*this->v[1] - this->r[1]*this->v[0]);
    double lsq = l_x*l_x + l_y*l_y + l_z*l_z;

    double r      = pow(this->distance(planet_2), 3.0);
    double rsq    = pow(this->distance(planet_2), 2.0);
    double mratio = planet_2.mass/consts::mass_sun;
    double delr   = this->r[dim] - planet_2.r[dim];
    double accN   = -(consts::four_pi_sq*mratio)*(delr/r);
    double acc    = accN*(1.0 + (3.0*lsq)/(rsq*consts::csq));
    return acc;
}

double planet::kinetic_energy() const
{
    double KE     = 0.5*this->mass*(pow(this->v[0], 2.0) + pow(this->v[1], 2.0) + pow(this->v[2], 2.0));
    return KE;
}

double planet::potential_energy(planet planet_2) const
{
    double dinom = (this->beta-1)*pow(this->distance(planet_2),this->beta-1);
    double PE    = - consts::G*this->mass*planet_2.mass/dinom;
    return PE;
}

planet::~planet()
{

}
