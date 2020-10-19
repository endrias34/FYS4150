# include <iostream>
# include <cmath>
# include "solver.h"

solver:: solver()
{
    planet* bodies;
    this->num_bodies = 0;
    this->fixed_sun  = true;
    this->bodies     = new planet[this->num_bodies];
    this->ofiles     = new std::ofstream[this->num_bodies];
    this->bodies     = bodies;

}

solver::solver(planet* bodies, int n, bool implicit_sun)
{
    this->num_bodies = n;
    this->fixed_sun  = implicit_sun;
    this->bodies     = new planet[num_bodies];
    this->ofiles     = new std::ofstream[num_bodies];

    for (int i = 0; i < num_bodies; i++)
    {
        this->bodies[i] = bodies[i];
    }
}

void solver::solve(double h, double t_max, double t_write, std::string method)
{
    int frame_write  = (int)(t_write/h + 0.5);       
    int total_frames = (int)(t_max/t_write + 0.5);
    int hmm          = abs(log10(h));
    int tnn          = t_max;
    std::string hstr = std::to_string(hmm);
    std::string tmm  = std::to_string(tnn);

    if (total_frames > 10000000){
        std::cout << "\nError: Too many lines writen out" << std::endl;
        std::cout << "Terminating program...\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->num_bodies; i++)
    {
        std::string filename = "/Users/endriasa/Projects/FYS4150/src/Project3/ObjectOrient/Results/";
        //filename.append("EscapeVr3_99_");
        filename.append("GR_");
        //filename.append("ThreeB_Msun_");
        filename.append(this->bodies[i].name);
        filename.append("_");
        filename.append(method);
        filename.append("_");
        filename.append(tmm);
        filename.append("_");
        filename.append(hstr);
        filename.append(".txt");

        this->ofiles[i].open(filename.c_str());
        this->ofiles[i] << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        this->ofiles[i] << "t x y z" << std::endl;   

        this->write_row_to_file(i, 0.0, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
    }

    /* ========================== Executing main algorithms ========================== */
   
    clock_t t_0 = clock(); 

    if (method.compare("euler") == 0)
    {
        this->euler(h, t_max, frame_write);
    }

    else if (method.compare("verlet") == 0)
    {
        this->verlet(h, t_max, frame_write);
    }

    else if(method.compare("verlet_GR") == 0)
    {
        this->verlet_GR(h, t_max, frame_write);
    }

    else {
        std::cout << "\nError: Invalid method! Choose either 'euler' or 'verlet' or 'verlet_GR'" << std::endl;
        std::cout << "Terminating program.." << std::endl;
        exit(EXIT_FAILURE);
    }

    clock_t t_1 = clock();      // Done timing
   double time_used = (double)(t_1 - t_0)/CLOCKS_PER_SEC;  
   
    /* Print time consumption */
    //std::cout << time_used << std::endl;

    for (int i = 0; i < this->num_bodies; i++){
        this->ofiles[i].close();   
    }
}

void solver::euler(double h, double t_max, int frame_write)
{
    planet* bodies_curr = new planet[this->num_bodies];
    int frame           = 0;
    double t            = 0.0;

    while (t <= t_max)
    {
        t += h;     
        frame++;    

       //this->compute_energy(t,frame_write, frame);
       //this->compute_angular_momentum(t,frame_write, frame);

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Copy of all planets for current time step
        }

        /* Solve differential equations for all bodies. */
        for (int i = 0; i < this->num_bodies; i++){

            /* Execute euler time-stepping for both/all directions. */
            for (int dim = 0; dim < cnst::num_dims; dim++){
                bodies_curr[i].a[dim] = this->compute_total_acc(bodies_curr[i], bodies_curr, dim);
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim];
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h*bodies_curr[i].a[dim];
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
            }
        }
    }

    delete[] bodies_curr;
}

/* Function that solves the n-body problem using the "Velocity Verlet" algorithm for time-stepping
 * the position and velocity in the x- and y-direction (can be extended to include z as well). The
 * only assumed knowlede about the body objects used in this function, is that they have a member
 * function compute_total_acc() that computes and returns the total acceleration on a body. */
void solver::verlet(double h, double t_max, int frame_write){
    planet* bodies_curr = new planet[this->num_bodies];
    int frame = 0;                      // Counting time steps for writing to file
    double h_squared_half = h*h/2.0;
    double h_half = h/2.0;
    double t = 0.0;

    /* Make copy of all bodies to track current and next time step. Also compute init. acceleration. */
    for (int i = 0; i < this->num_bodies; i++){
        bodies_curr[i] = this->bodies[i];

        for (int dim = 0; dim < cnst::num_dims; dim++){
            bodies_curr[i].a[dim] = this->compute_total_acc(this->bodies[i], this->bodies, dim);
        }
    }

    /* Loop until maximum times is reached. */
    while (t <= t_max)
    {
        t += h;     // Increase time by step size
        frame++;    // Count to next frame

        //this->compute_energy(t,frame_write, frame);
        //this->compute_angular_momentum(t,frame_write, frame);

        /* Update positions for all bodies in all directions. */
        for (int i = 0; i < this->num_bodies; i++){
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim] +
                        h_squared_half*bodies_curr[i].a[dim];

            }
        }

        /* New loop over bodies and directions to compute the updated values for acceleration
        * and then velocities. Can't do in same loop above since all r_{i+1} must be done for a_{i+1}. */
        for (int i = 0; i < this->num_bodies; i++){
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].a[dim] = this->compute_total_acc(this->bodies[i], this->bodies, dim);
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h_half*
                        (this->bodies[i].a[dim] + bodies_curr[i].a[dim]);
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
            }
        }

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Update current bodies
        }
    }
    delete[] bodies_curr;
}

/* Specific Verlet only for use when running mercury perihilion GR case. */
void solver::verlet_GR(double h, double t_max, int frame_write){
    planet* bodies_curr = new planet[this->num_bodies];
    int frame = 0;                      // Counting time steps for writing to file
    double h_squared_half = h*h/2.0;
    double h_half = h/2.0;
    double t = 0.0;
    double r_p[2];  // To hold r_perhilion for both GR mercury and ordinary mercury
    double x_p[2];  // To hold x_perhilion for both GR mercury and ordinary mercury
    double y_p[2];  // To hold y_perhilion for both GR mercury and ordinary mercury

    r_p[0] = 1.0;   // Definitely outside mercury orbit
    r_p[1] = 1.0;   // Definitely outside mercury orbit

    /* Make copy of all bodies to track current and next time step. Also compute init. acceleration. */
    for (int i = 0; i < this->num_bodies; i++){
        bodies_curr[i] = this->bodies[i];

        for (int dim = 0; dim < cnst::num_dims; dim++){
            bodies_curr[i].a[dim] = this->compute_total_acc_GR(this->bodies[i], dim);
        }
    }

    /* Loop until maximum times is reached. */
    while (t <= t_max){
        t += h;     // Increase time by step size
        frame++;    // Count to next frame


        /* Update positions for all bodies in all directions. */
        for (int i = 0; i < this->num_bodies; i++){
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim] +
                        h_squared_half*bodies_curr[i].a[dim];

            }
        }

        /* New loop over bodies and directions to compute the updated values for acceleration
        * and then velocities. Can't do in same loop above since all r_{i+1} must be done for a_{i+1}. */
        for (int i = 0; i < this->num_bodies; i++){
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].a[dim] = this->compute_total_acc_GR(this->bodies[i], dim);
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h_half*
                        (this->bodies[i].a[dim] + bodies_curr[i].a[dim]);
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
            }
        }

        // Extra stuff for perihelion computation
        if (t > (t_max - cnst::mercury_year)){
            for (int i = 0; i< this->num_bodies; i++){
               double x = this->bodies[i].r[0];
               double y = this->bodies[i].r[1];
               double z = this->bodies[i].r[2];
               double r = sqrt(x*x+y*y+z*z);

               /* Find minimum distance to sun. */
               if (r < r_p[i]){
                   r_p[i] = r;
                   x_p[i] = x;
                   y_p[i] = y;
               }
            }
        }

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Update current bodies
        }
    }

    for (int i =0; i < this->num_bodies; i++)
    {
        double theta_p = atan(y_p[i]/x_p[i])*(180.0/acos(-1))*3600.0;
        std::cout << "\n" << this->bodies[i].name << ": theta_p = " << theta_p
                  << ",  x_p = " << x_p[i] << ", y_p = " << y_p[i] << std::endl;
    }

    delete[] bodies_curr;
}


/* Function that, through calls to planet::compute_accel(), computes total acceleration of subject
 * (in dim-direction) due to the sum of gravitational forces from all planets in array objects. */
double solver::compute_total_acc(planet subject, planet* objects, int dim) const {
    double total_accel = 0.0;

    /* Simulation is being run with a fixed sun at the origin. */
    if (this->fixed_sun == true){
        planet sun("sun", cnst::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        total_accel += subject.compute_acc(sun, dim);     // Include force from sun separately
    }

    for (int i = 0; i < this->num_bodies; i++){
        if (objects[i].name.compare(subject.name) != 0){    // No force on itself
            //if (dim == 0) std::cout << "While using force from " << planets[i].name << std::endl;
            total_accel += subject.compute_acc(objects[i], dim);
        }
    }

    return total_accel;
}

/* Separate acceleration function for GR case. */
double solver::compute_total_acc_GR(planet subject, int dim) {
    double total_accel = 0.0;

    /* GR Simulation is being run with a fixed sun at the origin. */
    if (this->fixed_sun == true){
        planet sun("sun", cnst::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        if (subject.name.compare("mercury_GR") == 0){
            total_accel += subject.compute_acc_mercury_GR(sun, dim);    // Force for GR Mercury
        }

        else {
            total_accel += subject.compute_acc(sun, dim);   // Fore regular non-GR Mercury
        }
    }

    else {
        std::cout << "Error: Sun shall be fixed for the GR simulation!" << std::endl;
        std::cout << "Terminating program.." << std::endl;
        exit(EXIT_FAILURE);
    }


    return total_accel;
}

void solver::compute_energy(double time, int frame_write, int frame) const
{
    planet sun("sun", cnst::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    double K_E = 0.0;
    double P_E = 0.0;

    for (int i = 0; i < this->num_bodies; i++){
        K_E += this->bodies[i].compute_kinetic_energy();
        if (this->fixed_sun == true){
            P_E += this->bodies[i].compute_potential_energy(sun);
        }
        for (int j = 0; j < this->num_bodies; j++){
            if (j > i){
                P_E += this->bodies[i].compute_potential_energy(this->bodies[j]);
            }
        }
    }

    if (frame % frame_write == 0)
    {
     std::cout << std::setprecision(30) <<time <<" "<< K_E << " " << P_E << " " << K_E + P_E << std::endl;

    }

    
}


void solver::compute_angular_momentum(double time, int frame_write, int frame) const 
{
    double m, x, y, z, v_x, v_y, v_z;
    double l_x = 0.0, l_y = 0.0, l_z = 0.0;

    for (int i = 0; i < this->num_bodies; i++)
    {
        x = this->bodies[i].r[0];
        y = this->bodies[i].r[1];
        z = this->bodies[i].r[2];
        v_x = this->bodies[i].v[0];
        v_y = this->bodies[i].v[1];
        v_z = this->bodies[i].v[2];
        m = this->bodies[i].mass;
        l_x += m*(y*v_z - z*v_y);   // x-component of cross product
        l_y += m*(z*v_x - x*v_z);   // y-component of cross product
        l_z += m*(x*v_y - y*v_x);   // z-component of cross product
    }


    if (frame % frame_write == 0)
    {
     std::cout << std::setprecision(30) <<time <<" "<< l_x << " " << l_y << " " << l_z << " "<< sqrt(l_x*l_x+l_y*l_y+l_z*l_z) << std::endl;

    }
    
}

/* Function that displays the center of mass of the system. */
void solver::compute_center_mass() const {
    double M = 0.0;
    double r_x = 0.0, r_y = 0.0, r_z = 0.0;

    for (int i = 0; i < this->num_bodies; i++){
        M += this->bodies[i].mass;
        r_x += this->bodies[i].mass*this->bodies[i].r[0];
        r_y += this->bodies[i].mass*this->bodies[i].r[1];
        r_z += this->bodies[i].mass*this->bodies[i].r[2];
    }

    r_x = r_x/M;
    r_y = r_y/M;
    r_z = r_z/M;

    std::cout << "Center mass x: " << std::setprecision(8) << r_x << std::endl;
    std::cout << "Center mass y: " << std::setprecision(8) << r_y << std::endl;
    std::cout << "Center mass z: " << std::setprecision(8) << r_z << std::endl;
}

/* Function that writes a row of values (t, x, y) to data member output file
 * ofiles[file_index] where file_index is unique per body i the simulation. */
void solver::write_row_to_file(int file_index, double t, double x, double y, double z){
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << t;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << x;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << y;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << z << std::endl;
}

/* Destructor that deallocates memory from dynamically allocated member variables. */
solver::~solver(){
    delete[] this->bodies;
    delete[] this->ofiles;
}
