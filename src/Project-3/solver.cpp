# include <iostream>
# include <cmath>
# include "solver.h"

solver:: solver()
{
    planet* n_planets;
    this->num_planets    = 0;
    this->stationary_sun = true;
    this->n_planets      = new planet[this->num_planets];
    this->ofiles         = new std::ofstream[this->num_planets];
    this->n_planets      = n_planets;

}

solver::solver(planet* n_planets, int n, bool state_sun)
{
    this->num_planets    = n;
    this->stationary_sun = state_sun;
    this->n_planets      = new planet[num_planets];
    this->ofiles         = new std::ofstream[num_planets];

    for (int i = 0; i < num_planets; i++)
    {
        this->n_planets[i] = n_planets[i];
    }
}

void solver::simulate(double h, double t_max, double t_write, std::string method, std::string casenr)
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

    for (int i = 0; i < this->num_planets; i++)
    {
    std::string filename = "/Users/endriasa/Projects/FYS4150/src/Project3/Codes/Results/";
    std::string bmm      = std::to_string(this->n_planets[i].beta);
    filename.append(this->n_planets[i].name);
    filename.append("_");
    filename.append(method);
    filename.append("_");
    filename.append(casenr);
    filename.append("_");
    filename.append(bmm);
    filename.append("_");
    filename.append(tmm);
    filename.append("_");
    filename.append(hstr); 
    filename.append(".txt");

    this->ofiles[i].open(filename.c_str());
    this->ofiles[i] << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    this->ofiles[i] << "Time   X      Y       Z" << std::endl; 
    this->write_to_file(i, 0.0, this->n_planets[i].r[0], this->n_planets[i].r[1], this->n_planets[i].r[2]);

    }

    /* ========================== Executing main algorithms ========================== */
   
    if (method.compare("euler") == 0)
    {
        this->euler(h, t_max, frame_write);
    }
    else if (method.compare("eulerTime") == 0)
    {
        this->eulerTime(h, t_max);
    }
    else if (method.compare("eulerEnergy") == 0)
    {
        this->eulerEnergy(h, t_max, frame_write);
    }
    else if (method.compare("eulerAngmom") == 0)
    {
        this->eulerEnergy(h, t_max, frame_write);
    }
    else if (method.compare("verlet") == 0)
    {
        this->verlet(h, t_max, frame_write);
    }
    else if (method.compare("verletTime") == 0)
    {
        this->verletTime(h, t_max);
    }
    else if (method.compare("verletEnergy") == 0)
    {
        this->verletEnergy(h, t_max, frame_write);
    }
    else if (method.compare("verletAngmom") == 0)
    {
        this->verletAngmom(h, t_max, frame_write);
    }
    else if(method.compare("verlet_GR") == 0)
    {
        this->verlet_GR(h, t_max, frame_write);
    }

    for (int i = 0; i < this->num_planets; i++)
    {
        this->ofiles[i].close();  
    }

}

void solver::euler(double h, double t_max, int frame_write)
{
    
    planet* n_planets_curr = new planet[this->num_planets];
    int frame   = 0;
    double t    = 0.0;

    while (t <= t_max)
    {
        t     += h;     
        frame ++;    

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                n_planets_curr[i].a[dim]  = this->total_acc(n_planets_curr[i], n_planets_curr, dim);
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim];
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h*n_planets_curr[i].a[dim];
            }

            if (frame % frame_write == 0)
            {
                this->write_to_file(i, t, this->n_planets[i].r[0], this->n_planets[i].r[1], this->n_planets[i].r[2]);
            }
        }
    }

    delete[] n_planets_curr;
}

void solver::eulerTime(double h, double t_max)
{   
    clock_t t_0 = clock(); 
    planet* n_planets_curr = new planet[this->num_planets];
    double t    = 0.0;

    while (t <= t_max)
    {
        t += h;     

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                n_planets_curr[i].a[dim]  = this->total_acc(n_planets_curr[i], n_planets_curr, dim);
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim];
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h*n_planets_curr[i].a[dim];
            }
        }
    }

    clock_t t_1 = clock();      
    double time_used = (double)(t_1 - t_0)/CLOCKS_PER_SEC;  
    std::cout << time_used << std::endl;

    delete[] n_planets_curr;
}

void solver::eulerEnergy(double h, double t_max, int frame_write)
{
    planet* n_planets_curr = new planet[this->num_planets];
    int frame   = 0;
    double t    = 0.0;

    while (t <= t_max)
    {
        t     += h;     
        frame ++;    

       this->energy(t,frame_write, frame);

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                n_planets_curr[i].a[dim]  = this->total_acc(n_planets_curr[i], n_planets_curr, dim);
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim];
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h*n_planets_curr[i].a[dim];
            }

        }

    }

    delete[] n_planets_curr;
}
void solver::eulerAngmom(double h, double t_max, int frame_write)
{
    planet* n_planets_curr = new planet[this->num_planets];
    int frame   = 0;
    double t    = 0.0;

    while (t <= t_max)
    {
        t     += h;     
        frame ++;    

       this->angular_momentum(t,frame_write, frame);

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                n_planets_curr[i].a[dim]  = this->total_acc(n_planets_curr[i], n_planets_curr, dim);
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim];
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h*n_planets_curr[i].a[dim];
            }

        }

    }

    delete[] n_planets_curr;
}

void solver::verlet(double h, double t_max, int frame_write)
{
    planet* n_planets_curr = new planet[this->num_planets];
    int frame              = 0;                     
    double h_squared_half  = h*h/2.0;
    double h_half          = h/2.0;
    double t               = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        n_planets_curr[i] = this->n_planets[i];

        for (int dim = 0; dim < consts::num_dims; dim++)
        {
            n_planets_curr[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
        }
    }

    while (t <= t_max)
    {
        t += h;     
        frame++;    

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim] +
                        h_squared_half*n_planets_curr[i].a[dim];

            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h_half*
                        (this->n_planets[i].a[dim] + n_planets_curr[i].a[dim]);
            }

            if (frame % frame_write == 0)
            {
                this->write_to_file(i, t, this->n_planets[i].r[0], this->n_planets[i].r[1], this->n_planets[i].r[2]);
            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }
    }

    delete[] n_planets_curr;
}

void solver::verletTime(double h, double t_max)
{
    clock_t t_0 = clock();
    planet* n_planets_curr = new planet[this->num_planets];                    
    double h_squared_half  = h*h/2.0;
    double h_half          = h/2.0;
    double t               = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        n_planets_curr[i] = this->n_planets[i];

        for (int dim = 0; dim < consts::num_dims; dim++)
        {
            n_planets_curr[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
        }
    }

    while (t <= t_max)
    {
        t += h;     

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim] +
                        h_squared_half*n_planets_curr[i].a[dim];

            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h_half*
                        (this->n_planets[i].a[dim] + n_planets_curr[i].a[dim]);
            }

        }

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }
    }

    clock_t t_1 = clock();      
    double time_used = (double)(t_1 - t_0)/CLOCKS_PER_SEC;  
    std::cout << time_used << std::endl;

    delete[] n_planets_curr;
}

void solver::verletEnergy(double h, double t_max, int frame_write)
{
    planet* n_planets_curr = new planet[this->num_planets];
    int frame              = 0;                     
    double h_squared_half  = h*h/2.0;
    double h_half          = h/2.0;
    double t               = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        n_planets_curr[i] = this->n_planets[i];

        for (int dim = 0; dim < consts::num_dims; dim++)
        {
            n_planets_curr[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
        }
    }

    while (t <= t_max)
    {
        t += h;     
        frame++;    

        this->energy(t,frame_write, frame);

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim] +
                        h_squared_half*n_planets_curr[i].a[dim];

            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h_half*
                        (this->n_planets[i].a[dim] + n_planets_curr[i].a[dim]);
            }

        }

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }
    }

    delete[] n_planets_curr;
}
void solver::verletAngmom(double h, double t_max, int frame_write)
{
    planet* n_planets_curr = new planet[this->num_planets];
    int frame              = 0;                     
    double h_squared_half  = h*h/2.0;
    double h_half          = h/2.0;
    double t               = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        n_planets_curr[i] = this->n_planets[i];

        for (int dim = 0; dim < consts::num_dims; dim++)
        {
            n_planets_curr[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
        }
    }

    while (t <= t_max)
    {
        t += h;     
        frame++;    

        this->angular_momentum(t,frame_write, frame);

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim] +
                        h_squared_half*n_planets_curr[i].a[dim];

            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].a[dim] = this->total_acc(this->n_planets[i], this->n_planets, dim);
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h_half*
                        (this->n_planets[i].a[dim] + n_planets_curr[i].a[dim]);
            }

        }

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];       
        }
    }

    delete[] n_planets_curr;
}

void solver::verlet_GR(double h, double t_max, int frame_write)
{
    planet* n_planets_curr = new planet[this->num_planets];
    int frame              = 0;                      
    double h_squared_half  = h*h/2.0;
    double h_half          = h/2.0;
    double t               = 0.0;

    double r_p[2];  
    double x_p[2];  
    double y_p[2];  

    r_p[0] = 1.0;   
    r_p[1] = 1.0;  

    for (int i = 0; i < this->num_planets; i++)
    {
        n_planets_curr[i] = this->n_planets[i];

        for (int dim = 0; dim < consts::num_dims; dim++)
        {
            n_planets_curr[i].a[dim] = this->total_acc_GR(this->n_planets[i], dim);
        }
    }

    while (t <= t_max)
    {
        t    += h;     
        frame++;   

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].r[dim] = n_planets_curr[i].r[dim] + h*n_planets_curr[i].v[dim] +
                        h_squared_half*n_planets_curr[i].a[dim];
            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            for (int dim = 0; dim < consts::num_dims; dim++)
            {
                this->n_planets[i].a[dim] = this->total_acc_GR(this->n_planets[i], dim);
                this->n_planets[i].v[dim] = n_planets_curr[i].v[dim] + h_half*
                        (this->n_planets[i].a[dim] + n_planets_curr[i].a[dim]);
            }

            if (frame % frame_write == 0)
            {
                this->write_to_file(i, t, this->n_planets[i].r[0], this->n_planets[i].r[1], this->n_planets[i].r[2]);
            }
        }

        if (t > (t_max - consts::mercury_year))
        {
            for (int i = 0; i< this->num_planets; i++)
            {
               double x = this->n_planets[i].r[0];
               double y = this->n_planets[i].r[1];
               double z = this->n_planets[i].r[2];
               double r = sqrt(x*x+y*y+z*z);

               if (r < r_p[i])
               {
                   r_p[i] = r;
                   x_p[i] = x;
                   y_p[i] = y;
               }

                double theta_p = atan(y_p[i]/x_p[i])*(180.0/consts::pi)*3600.0;
                std::cout  << t << " " << theta_p << " " << x_p[i] << " " << y_p[i]  << std::endl;

            }
        }

        for (int i = 0; i < this->num_planets; i++)
        {
            n_planets_curr[i] = this->n_planets[i];        
        }
    }

    delete[] n_planets_curr;
}

double solver::total_acc(planet subject, planet* objects, int dim) const
{
    double total_accel = 0.0;

    if (this->stationary_sun == true)
    {
        planet sun("sun", consts::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,2.0);
        total_accel += subject.acceleration(sun, dim);     
    }

    for (int i = 0; i < this->num_planets; i++)
    {
        if (objects[i].name.compare(subject.name) != 0)
        {    
            total_accel += subject.acceleration(objects[i], dim);
        }
    }

    return total_accel;
}

double solver::total_acc_GR(planet subject, int dim) 
{
    double total_accel = 0.0;

    if (this->stationary_sun == true)
    {
        planet sun("sun", consts::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0);

        if (subject.name.compare("mercury_GR") == 0)
        {
            total_accel += subject.acc_mercury_GR(sun, dim);   
        }

        else 
        {
            total_accel += subject.acceleration(sun, dim);   
        }
    }

    else {
        std::cout << "Error: Sun shall be stationary for the GR simulation!" << std::endl;
        std::cout << "Terminating program.." << std::endl;
        exit(EXIT_FAILURE);
    }


    return total_accel;
}

void solver::energy(double time, int frame_write, int frame) const
{
    planet sun("sun", consts::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,2.0);
    double K_E = 0.0;
    double P_E = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        K_E += this->n_planets[i].kinetic_energy();

        if (this->stationary_sun == true)
        {
            P_E += this->n_planets[i].potential_energy(sun);
        }

        for (int j = 0; j < this->num_planets; j++)
        {
            if (j > i)
            {
                P_E += this->n_planets[i].potential_energy(this->n_planets[j]);
            }
        }
    }

    if (frame % frame_write == 0)
    {
     std::cout << std::setprecision(30) <<time <<" "<< K_E << " " << P_E << " " << K_E + P_E << std::endl;

    }
    
}


void solver::angular_momentum(double time, int frame_write, int frame) const
{
    double m, x, y, z, v_x, v_y, v_z;
    double l_x = 0.0, l_y = 0.0, l_z = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        x = this->n_planets[i].r[0];
        y = this->n_planets[i].r[1];
        z = this->n_planets[i].r[2];

        v_x = this->n_planets[i].v[0];
        v_y = this->n_planets[i].v[1];
        v_z = this->n_planets[i].v[2];

        m = this->n_planets[i].mass;

        l_x += m*(y*v_z - z*v_y);   
        l_y += m*(z*v_x - x*v_z);   
        l_z += m*(x*v_y - y*v_x);   
    }

    if (frame % frame_write == 0)
    {
     std::cout << std::setprecision(30) <<time <<" "<< l_x << " " << l_y << " " << l_z << " "<< sqrt(l_x*l_x+l_y*l_y+l_z*l_z) << std::endl;

    }
    
}

void solver::center_of_mass() const
{
    double M   = 0.0;
    double r_x = 0.0, r_y = 0.0, r_z = 0.0;

    for (int i = 0; i < this->num_planets; i++)
    {
        M   += this->n_planets[i].mass;
        r_x += this->n_planets[i].mass*this->n_planets[i].r[0];
        r_y += this->n_planets[i].mass*this->n_planets[i].r[1];
        r_z += this->n_planets[i].mass*this->n_planets[i].r[2];
    }

    r_x = r_x/M;
    r_y = r_y/M;
    r_z = r_z/M;

    std::cout << "Center of mass x: " << std::setprecision(8) << r_x << std::endl;
    std::cout << "Center of mass y: " << std::setprecision(8) << r_y << std::endl;
    std::cout << "Center of mass z: " << std::setprecision(8) << r_z << std::endl;
}

void solver::write_to_file(int file_index, double t, double x, double y, double z)
{
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << t;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << x;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << y;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << z << std::endl;
}

solver::~solver()
{
    delete[] this->n_planets;
    delete[] this->ofiles;

}
