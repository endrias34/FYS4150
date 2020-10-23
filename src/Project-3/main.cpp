# include <iostream>
# include <cmath>
# include "planet.h"
# include "solver.h"
# include "All_const.h"

/* Commande line input structure:
* [code.exe] [case_sel] [method] [max. time] [Step size] [tw_fact] [beta]
* [case_sel]    --> case1(2body), case2(3body), case3(3body moving sun), case4(full solar system), case5(GR_mercury)
* [method]      --> verlet, euler 
* [max. time]   --> 100 (in years)
* [Step size]   --> 10^(-[step size]) (in years)
* [tw_fact]     --> factor for writing every tw_fact*h result
* [beta]        --> factor for Gravitational force
* */

void sun_initial_velocity(planet*, double, double&, double&, double&);

int main(int argc, char const *argv[])
{
    // define input parameters
    std::string case_sel = argv[1];  
    std::string method   = argv[2];
    double t_max         = atof(argv[3]);
    double h             = pow(10, -atof(argv[4]));
    int tw_fact          = atoi(argv[5]);
    double t_write       = tw_fact*h;
    double beta          = atof(argv[6]);

    planet* n_planets;
    bool stationary_sun;
    int num_planets;
        
    // Circular orbit 
    if (case_sel.compare("case1") == 0)
    {
        num_planets     = 1;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2.0*consts::pi, 0.0, beta);
        n_planets[0] = earth;

        solver earth_sun(n_planets, num_planets, stationary_sun);      
        earth_sun.simulate(h, t_max, t_write, method, case_sel); 
    }

    // Elliptical orbit with modified gravity
    else if (case_sel.compare("case11") == 0)
    {
        num_planets     = 1;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 5, 0.0, beta);
        n_planets[0] = earth;

        solver earth_sun(n_planets, num_planets, stationary_sun);      
        earth_sun.simulate(h, t_max, t_write, method, case_sel); 
    }

    // Escape velocity
    else if (case_sel.compare("case111") == 0)
    {
        num_planets     = 1;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2.0*consts::pi*sqrt(2), 0.0, beta);
        n_planets[0] = earth;

        solver earth_sun(n_planets, num_planets, stationary_sun);      
        earth_sun.simulate(h, t_max, t_write, method, case_sel); 

    }
    // Escape velocity*99%
    else if (case_sel.compare("case111b") == 0)
    {
        num_planets     = 1;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2.0*consts::pi*sqrt(2)*0.99, 0.0, beta);
        n_planets[0] = earth;

        solver earth_sun(n_planets, num_planets, stationary_sun);      
        earth_sun.simulate(h, t_max, t_write, method, case_sel); 

    }
        // Escape velocity for 1/r^3
    else if (case_sel.compare("case111c") == 0)
    {
        num_planets     = 1;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2.0*consts::pi*1.0001, 0.0, beta);
        n_planets[0] = earth;

        solver earth_sun(n_planets, num_planets, stationary_sun);      
        earth_sun.simulate(h, t_max, t_write, method, case_sel); 

    }
        // Escape velocity*99% of 1/r^3
    else if (case_sel.compare("case111d") == 0)
    {
        num_planets     = 1;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2.0*consts::pi*sqrt(2), 0.0, beta);
        n_planets[0] = earth;

        solver earth_sun(n_planets, num_planets, stationary_sun);      
        earth_sun.simulate(h, t_max, t_write, method, case_sel); 

    }
    else if (case_sel.compare("case3") == 0)
    {
        num_planets     = 2;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 9.607870614931846E-01,  2.732310489311487E-01, -1.615338753099692E-05,
        -4.979870116841902E-03*365,  1.648427775007217E-02*365, -1.352988249682294E-06*365, beta);
        planet jupiter("jupiter", 1.9E+27, 2.511420852264080E+00, -4.457270905787417E+00, -3.769750456705974E-02,  
        6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        n_planets[0] = earth;
        n_planets[1] = jupiter;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
    }
    else if (case_sel.compare("case31") == 0)
    {
        num_planets     = 2;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 9.607870614931846E-01,  2.732310489311487E-01, -1.615338753099692E-05,
        -4.979870116841902E-03*365,  1.648427775007217E-02*365, -1.352988249682294E-06*365, beta);
        planet jupiter("jupiter", 10*1.9E+27, 2.511420852264080E+00, -4.457270905787417E+00, -3.769750456705974E-02,  
        6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        n_planets[0] = earth;
        n_planets[1] = jupiter;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
    }
    else if (case_sel.compare("case311") == 0)
    {
        num_planets     = 2;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet earth("earth", 6.0E+24, 9.607870614931846E-01,  2.732310489311487E-01, -1.615338753099692E-05,
        -4.979870116841902E-03*365,  1.648427775007217E-02*365, -1.352988249682294E-06*365, beta);
        planet jupiter("jupiter", 1000*1.9E+27, 2.511420852264080E+00, -4.457270905787417E+00, -3.769750456705974E-02,  
        6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        n_planets[0] = earth;
        n_planets[1] = jupiter;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
    }

    else if (case_sel.compare("case4") == 0)
    {
        num_planets     = 3;
        stationary_sun  = false;
        n_planets       = new planet[num_planets];

        double cmx = -0.0036645731;
        double cmy = 0.0022203996;
        double cmz = 5.1546204e-05;

        planet earth("earth", 6.0E+24, 9.607870614931846E-01-cmx,  2.732310489311487E-01-cmy, -1.615338753099692E-05-cmz,
        -4.979870116841902E-03*365,  1.648427775007217E-02*365, -1.352988249682294E-06*365, beta);

        planet jupiter("jupiter", 1.9E+27, 2.511420852264080E+00-cmx, -4.457270905787417E+00-cmy, -3.769750456705974E-02-cmz,  
        6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        n_planets[0]  = earth;
        n_planets[1]  = jupiter;

        double v_x, v_y, v_z;
        sun_initial_velocity(n_planets, num_planets, v_x, v_y, v_z);
        planet sun("sun", 2.0E+30, -6.056797627705779E-03-cmx,  6.456103347534338E-03-cmy,  8.740800510641314E-05-cmz,
         v_x, v_y, v_z, beta);
        n_planets[2]  = sun;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
        //earth_jupiter_sun.center_of_mass();

    }

    else if (case_sel.compare("case41") == 0)
    {
        num_planets     = 3;
        stationary_sun  = false;
        n_planets       = new planet[num_planets];

        double cmx = 0.017636979;
        double cmy = -0.035549325;
        double cmz = -0.00026816992;

        planet earth("earth", 6.0E+24, 9.607870614931846E-01-cmx,  2.732310489311487E-01-cmy, -1.615338753099692E-05-cmz,
        -4.979870116841902E-03*365,  1.648427775007217E-02*365, -1.352988249682294E-06*365, beta);

        planet jupiter("jupiter", 10*1.9E+27, 2.511420852264080E+00-cmx, -4.457270905787417E+00-cmy, -3.769750456705974E-02-cmz,  
        6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        n_planets[0]  = earth;
        n_planets[1]  = jupiter;

        double v_x, v_y, v_z;
        sun_initial_velocity(n_planets, num_planets, v_x, v_y, v_z);
        planet sun("sun", 2.0E+30, -6.056797627705779E-03-cmx,  6.456103347534338E-03-cmy,  8.740800510641314E-05-cmz,
         v_x, v_y, v_z, beta);
        n_planets[2]  = sun;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
        //earth_jupiter_sun.center_of_mass();

    }

    else if (case_sel.compare("case411") == 0)
    {
        num_planets     = 3;
        stationary_sun  = false;
        n_planets       = new planet[num_planets];

        double cmx = 1.2204063;
        double cmy = -2.1681764;
        double cmz = -0.018320598;

        planet earth("earth", 6.0E+24, 9.607870614931846E-01-cmx,  2.732310489311487E-01-cmy, -1.615338753099692E-05-cmz,
        -4.979870116841902E-03*365,  1.648427775007217E-02*365, -1.352988249682294E-06*365, beta);

        planet jupiter("jupiter", 1000*1.9E+27, 2.511420852264080E+00-cmx, -4.457270905787417E+00-cmy, -3.769750456705974E-02-cmz,  
        6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        n_planets[0]  = earth;
        n_planets[1]  = jupiter;

        double v_x, v_y, v_z;
        sun_initial_velocity(n_planets, num_planets, v_x, v_y, v_z);
        planet sun("sun", 2.0E+30, -6.056797627705779E-03-cmx,  6.456103347534338E-03-cmy,  8.740800510641314E-05-cmz,
         v_x, v_y, v_z, beta);
        n_planets[2]  = sun;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
        //earth_jupiter_sun.center_of_mass();

    }

    else if (case_sel.compare("case5") == 0)
    {
        num_planets     = 10;
        stationary_sun  = false;
        n_planets       = new planet[num_planets];

        double cmx = -6.2806658e-05;
        double cmy = 0.0001188617;
        double cmz = 7.2648688e-07;

        planet mercury("mercury", 3.3E+23, 2.939793096020483E-01 -cmx, -2.743915003505926E-01 -cmy, -5.038494575001255E-02 -cmz,
          1.364847611039539E-02*365,  2.186100432554229E-02*365,  5.343438290447355E-04*365, beta);

        planet venus("venus", 4.9E+24, -8.198092990544088E-02 -cmx,  7.216213946591721E-01 -cmy,  1.428289213118692E-02 -cmz,
         -2.018960022294621E-02*365, -2.250708580201211E-03*365,  1.134043425643825E-03*365, beta);
        
        planet earth("earth", 6.0E+24, 9.547302638654790E-01 -cmx,  2.796871522786831E-01 -cmy,  7.125461757541623E-05 -cmz,
         -4.987197577527832E-03*365,  1.647924718197501E-02*365, -1.134205905479237E-06*365, beta);

        planet mars("mars", 6.6E+23, 1.338171629453269E+00 -cmx,  4.375797272123467E-01 -cmy, -2.385380064615869E-02 -cmz,
         -3.744983310791084E-03*365,  1.451508475315047E-02*365,  3.961953348278991E-04*365, beta);

        planet jupiter("jupiter", 1.9E+27, 2.511420852264080E+00 -cmx, -4.457270905787417E+00 -cmy, -3.769750456705974E-02 -cmz,
          6.481403418791627E-03*365,  4.061705749755052E-03*365, -1.618678241608630E-04*365, beta);

        planet saturn("saturn", 5.5E+26, 5.112689600325091E+00 -cmx, -8.583651246796396E+00 -cmy, -5.429222143028824E-02 -cmz,
          4.483016338695055E-03*365,  2.840107438540901E-03*365, -2.277544695614591E-04*365, beta);

        planet uranus("uranus", 8.8E+25, 1.555207444103421E+01 -cmx,  1.222492019985741E+01 -cmy, -1.560755393502417E-01 -cmz,
         -2.459615755280528E-03*365,  2.908916217372186E-03*365,  4.275096308560520E-05*365, beta);
        
        planet neptune("neptune", 1.03E+26, 2.940840887606790E+01 -cmx, -5.486811564721911E+00 -cmy, -5.647567541576598E-01 -cmz,
          5.547846555658923E-04*365,  3.104101465188291E-03*365, -7.699814259230016E-05*365, beta);

        planet pluto("pluto", 1.31E+22, 1.380705794972566E+01 -cmx, -3.120338950577318E+01 -cmy, -6.548628267867134E-01 -cmz,
          2.949228879130654E-03*365,  6.086689603034977E-04*365, -9.081761772829586E-04*365, beta);

        n_planets[0] = mercury;
        n_planets[1] = venus;
        n_planets[2] = earth;
        n_planets[3] = mars;
        n_planets[4] = jupiter;
        n_planets[5] = saturn;
        n_planets[6] = uranus;
        n_planets[7] = neptune;
        n_planets[8] = pluto;

        double v_x, v_y, v_z;
        sun_initial_velocity(n_planets, num_planets, v_x, v_y, v_z);
        planet sun("sun", 2.0E+30, -6.056797627705779E-03 -cmx,  6.456103347534338E-03 -cmy,  8.740800510641314E-05 -cmz,
         v_x, v_y, v_z, beta);
        n_planets[9]  = sun;

        solver earth_jupiter_sun(n_planets, num_planets, stationary_sun);      
        earth_jupiter_sun.simulate(h, t_max, t_write, method, case_sel); 
        //earth_jupiter_sun.center_of_mass();

    }
    
    else if (case_sel.compare("case6") == 0)
    {
        num_planets     = 2;
        stationary_sun  = true;
        n_planets       = new planet[num_planets];

        planet mercury_GR("mercury_GR", 3.3E23, 0.3075, 0.0 , 0.0, 0.0, 12.44, 0.0,beta);
        planet mercury("mercury", 3.3E23, 0.3075, 0.0, 0.0, 0.0, 12.44, 0.0,beta);

        n_planets[0] = mercury_GR;
        n_planets[1] = mercury;

        solver mercury_sun(n_planets, num_planets, stationary_sun);      
        mercury_sun.simulate(h, t_max, t_write, method, case_sel); 

    }

    delete[] n_planets;
    return 0;
}

void sun_initial_velocity(planet* n_planets, double num_planets, double& v_x, double& v_y, double& v_z)
{
    double tot_mom_x = 0.0;
    double tot_mom_y = 0.0;
    double tot_mom_z = 0.0;

    for (int i = 0; i < num_planets-1; i++){
        tot_mom_x += n_planets[i].mass*n_planets[i].v[0];
        tot_mom_y += n_planets[i].mass*n_planets[i].v[1];
        tot_mom_z += n_planets[i].mass*n_planets[i].v[2];
    }

    v_x = -tot_mom_x/consts::mass_sun;
    v_y = -tot_mom_y/consts::mass_sun;
    v_z = -tot_mom_z/consts::mass_sun;
}