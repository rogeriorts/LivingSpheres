#pragma once
#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>
#include "CImg.h"
#include "util.h"
#include <algorithm>
#include "grid.h"
#include <chrono>
#include <omp.h>


#define M_PI 3.14159265358979323846

using namespace cimg_library;


struct cell
{
    double x;
    double y;
    double vx;
    double vy;
    std::vector<double> fx;
    std::vector<double> fy;
    double k;
    double mass;
    double damp_ratio;
    int is_wall;
    std::vector<int> contacts;
    std::vector<double> contact_vector_x;
    std::vector<double> contact_vector_y;
    int number_of_contacts = 0;
    int grid_block = -1;

    void add_contact(int cell_b_, double contact_vector_x_,
                     double contact_vector_y_)
    {
        number_of_contacts++;
        
        if (number_of_contacts - 1 >= contacts.size())
        {
            contacts.push_back(cell_b_);
            contact_vector_x.push_back(contact_vector_x_);
            contact_vector_y.push_back(contact_vector_y_);
        }
        else
        {
            contacts[number_of_contacts - 1] = cell_b_;
            contact_vector_x[number_of_contacts - 1] = contact_vector_x_;
            contact_vector_y[number_of_contacts - 1] = contact_vector_y_;
        }
        
    }
};


struct cells_type
{

    std::vector<cell> cells_collection;
    int number_of_cells = 0;
    
    
};



bool verify_boundary_contact(double cell_position, double boundary_max, double radius, double &boundary_position)
{
    bool is_contact = false;

    if (cell_position < radius)
    {
        boundary_position = 0.0;
        is_contact = true;
    }
    else if (cell_position > boundary_max - radius)
    {
        boundary_position = boundary_max;
        is_contact = true;
    }
    return is_contact;
}

void add_cells_to_grid(cells_type &cells, grid_type &grid)
{
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        cells.cells_collection[i].grid_block = grid.add_cell_to_block(
            cells.cells_collection[i].x, cells.cells_collection[i].y, i);
        
    }
}

void find_contact_between_pair_of_cells(
    cells_type &cells, int i, int j, double radius)
{
    if ((cells.cells_collection[i].is_wall > -1 && cells.cells_collection[j].is_wall > -1) ||
         j <= i)
        return;

    double dist_x = cells.cells_collection[j].x - cells.cells_collection[i].x;
    double dist_y = cells.cells_collection[j].y - cells.cells_collection[i].y;
    double squared_distance = dist_x * dist_x + dist_y * dist_y;

    if (squared_distance < 4.0 * radius * radius)
    {
        cells.cells_collection[i].add_contact(j, dist_x, dist_y);
    }
}

void find_contacts_grid(cells_type &cells, grid_type &grid,
                        double radius, bool has_boundaries, double width = 500.0, double height = 500.0)
{
    for (size_t i = 0; i < cells.number_of_cells; i++)
    {
        cells.cells_collection[i].number_of_contacts = 0;
    }

    #pragma omp parallel for
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        int grid_k = cells.cells_collection[i].grid_block;
        int grid_j = grid_k / grid.nx; 
        int grid_i = grid_k - grid_j * grid.nx;

        for(int grid_i_neighbor = std::max(0,grid_i - 1); 
            grid_i_neighbor <= std::min(grid.nx - 1,grid_i + 1); 
            grid_i_neighbor++
        )
        {
            for(int grid_j_neighbor = std::max(0,grid_j - 1); 
                grid_j_neighbor <= std::min(grid.ny - 1,grid_j + 1); 
                grid_j_neighbor++
            )
            {
                int k_neighbor = grid_i_neighbor + grid.nx * grid_j_neighbor;
                for (int c_b = 0; c_b < grid.blocks[k_neighbor].next_position_on_block; c_b++)
                {
                    int j = grid.blocks[k_neighbor].cell_ids[c_b];
                    find_contact_between_pair_of_cells(cells, i, j, radius);
                }
            }
        }

        if (has_boundaries)
        {
            double boundary_position;
            // verify horizontal contacts with boundaries
            if (verify_boundary_contact(cells.cells_collection[i].x, width, radius, boundary_position))
            {
                double dist_x = boundary_position - cells.cells_collection[i].x;
                cells.cells_collection[i].add_contact(-1, dist_x, 0.0);                        
                
            }
            if (verify_boundary_contact(cells.cells_collection[i].y, height, radius, boundary_position))
            {
                double dist_y = boundary_position - cells.cells_collection[i].y;
                cells.cells_collection[i].add_contact(-1, 0.0, dist_y);
            }
        }
    }
    
        
    
}

void reset_forces(cells_type &cells)
{
    for(int i = 0 ; i< cells.number_of_cells ; i++)
    {
        std::fill(cells.cells_collection[i].fx.begin(), 
        cells.cells_collection[i].fx.end(), 0.0);

        std::fill(cells.cells_collection[i].fy.begin(), 
        cells.cells_collection[i].fy.end(), 0.0);
    }
}

void calculate_force_due_contacts(cells_type &cells, double radius, double friction_coefficient)
{
    #pragma omp parallel for
    for (int cell = 0; cell < cells.number_of_cells; cell++)
    {
        int t = omp_get_thread_num();

        for(int contact = 0;contact < cells.cells_collection[cell].number_of_contacts; contact++)
        {
            double dist_x = cells.cells_collection[cell].contact_vector_x[contact];
            double dist_y = cells.cells_collection[cell].contact_vector_y[contact];
            double squared_distance = dist_x * dist_x + dist_y * dist_y;

            // calculate normal vector between the contacting cells
            double distance = sqrt(squared_distance);
            double x_normal = dist_x / distance;
            double y_normal = dist_y / distance;

            double x_tangential = y_normal;
            double y_tangential = -x_normal;

            int i = cell;
            int j = cells.cells_collection[cell].contacts[contact];

            double avg_k, avg_mass, avg_damp_ratio, vx_b, vy_b, spring_deformation;

            if (j < 0)
            {
                avg_k = cells.cells_collection[i].k;
                avg_mass = cells.cells_collection[i].mass;
                avg_damp_ratio = cells.cells_collection[i].damp_ratio;
                vx_b = 0.0;
                vy_b = 0.0;
                spring_deformation = radius - distance;
            }
            else
            {
                avg_k = (cells.cells_collection[i].k + cells.cells_collection[j].k) * 0.5;
                avg_mass = (cells.cells_collection[i].mass + cells.cells_collection[j].mass) * 0.5;
                avg_damp_ratio = (cells.cells_collection[i].damp_ratio + cells.cells_collection[j].damp_ratio) * 0.5;
                vx_b = cells.cells_collection[j].vx;
                vy_b = cells.cells_collection[j].vy;
                spring_deformation = 2.0 * radius - distance;
            }
            // spring force
            if(spring_deformation < 0.0) continue;

            double spring_force = std::min(avg_k * spring_deformation, avg_k * radius);
            
            // damping force

            double vx_rel = vx_b - cells.cells_collection[i].vx;
            double vy_rel = vy_b - cells.cells_collection[i].vy;

            double damp_coeff = 2.0 * avg_damp_ratio * sqrt(avg_mass * radius * avg_k);

            double v_normal = vx_rel * x_normal + vy_rel * y_normal;
            double v_tangential = vx_rel * x_tangential + vy_rel * y_tangential;

            double vx_normal = v_normal * x_normal;
            double vy_normal = v_normal * y_normal;

            double vx_tangential = v_tangential * x_tangential;
            double vy_tangential = v_tangential * y_tangential;

            double damping_force_x = -vx_normal * damp_coeff;
            double damping_force_y = -vy_normal * damp_coeff;

            double abs_tangential = abs(v_tangential);
            double tangential_force_x = 0.0;
            double tangential_force_y = 0.0;
            if (abs_tangential > 1e-15)
            {
                tangential_force_x = -friction_coefficient * spring_force * vx_tangential / abs_tangential;
                tangential_force_y = -friction_coefficient * spring_force * vy_tangential / abs_tangential;
            }

            double normal_force_x = spring_force * x_normal;
            double normal_force_y = spring_force * y_normal;

            // calculate new velocities

            double fx = -normal_force_x - damping_force_x - tangential_force_x;
            double fy = -normal_force_y - damping_force_y - tangential_force_y;

            cells.cells_collection[i].fx[t] += fx;
            cells.cells_collection[i].fy[t] += fy;
            if (j >= 0)
            {
                cells.cells_collection[j].fx[t] -= fx;
                cells.cells_collection[j].fy[t] -= fy;
            }
        }
    }
}

void compute_new_position(cells_type &cells, double timestep, double gravity, int number_of_threads)
{
    #pragma omp parallel for
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        // gravity influence
        if (cells.cells_collection[i].is_wall < 0)
        {
            // calculating accelerations
            double fx=0.0;
            double fy=0.0;
            
            for (int t=0;t<number_of_threads ; t++)
            {
                fx += cells.cells_collection[i].fx[t];
                fy += cells.cells_collection[i].fy[t];
            }

            double ax = fx / cells.cells_collection[i].mass;
            double ay = fy / cells.cells_collection[i].mass + gravity;

            // calculate new positions
            //  s = s0 + v0 * t + a * t^2 /2
            double acceleration_timeste_factor = timestep * timestep * 0.5;
            cells.cells_collection[i].x += cells.cells_collection[i].vx * timestep + ax * acceleration_timeste_factor;
            cells.cells_collection[i].y += cells.cells_collection[i].vy * timestep + ay * acceleration_timeste_factor;

            // calculating new velocities
            // v = v0 + a * t
            cells.cells_collection[i].vx += ax * timestep;
            cells.cells_collection[i].vy += ay * timestep;
        }
    }
}

using namespace std::chrono;

bool search_trigger(cells_type &cells, double max_velocity_search)
{
    
    double squared_max_velocity = max_velocity_search*max_velocity_search;
    for(int i =0; i<cells.number_of_cells;i++)
    {   
        if (cells.cells_collection[i].vx*cells.cells_collection[i].vx + 
            cells.cells_collection[i].vy*cells.cells_collection[i].vy > squared_max_velocity)
            return true;
    }
    return false;
}

void do_one_iteration(cells_type &cells, grid_type &grid,
                      double radius, double timestep, double width, double height, double friction_coefficient, double gravity,
                      bool has_boundaries,
                      double &time_finding_contacts,
                      double &time_in_forces, double &time_in_move, int &n_timesteps, int &n_searchs,double max_velocity_search,
                      int number_of_threads)
{
    auto start = high_resolution_clock::now();
    
    
    if(search_trigger(cells, max_velocity_search))
    {
        
            grid.reset_grid();
            add_cells_to_grid(cells, grid);
            find_contacts_grid(cells, grid, radius, has_boundaries, width, height);
            n_searchs++;
        
    }
    auto stop = high_resolution_clock::now();
    time_finding_contacts += duration<double>(stop - start).count();

    start = high_resolution_clock::now();
    reset_forces(cells);
    calculate_force_due_contacts(cells, radius, friction_coefficient);
    stop = high_resolution_clock::now();
    time_in_forces += duration<double>(stop - start).count();
    
    start = high_resolution_clock::now();
    compute_new_position(cells, timestep, gravity,number_of_threads);
    stop = high_resolution_clock::now();
    time_in_move += duration<double>(stop - start).count();
    n_timesteps++;
}

struct Environment
{
    cells_type cells;
    grid_type grid;

    double radius;
    double friction_coefficient;
    double initial_number_of_cells;
    double width;
    double height;
    double gravity;
    bool cells_initialized;
    bool walls_initialized;

    std::vector<double> wall_ids;
    int number_of_walls;

    double has_boundaries;

    double timestep;

    double search_multiplier;

    Environment(
        double radius_,
        double width_,
        double height_,
        double grid_size_multiplier_,
        double search_multiplier_,
        double friction_coefficient_,
        bool has_boundaries_ = true,
        double gravity_ = -9.81)
    {
        radius = radius_;
        height = height_;
        width = width_;
        friction_coefficient = friction_coefficient_;
        search_multiplier = search_multiplier_;
        gravity = gravity_;
        has_boundaries = has_boundaries_;
        cells_initialized = false;
        walls_initialized = false;
        number_of_walls = 0;
        grid.grid_initialize(radius, width, height, grid_size_multiplier_);
    }

    void initialize_living_cells(int initial_number_of_cells,
                                 double k_min, double k_max,
                                 double damping_min, double damping_max, double x_min, double y_min,
                                 double x_max, double y_max, int number_of_threads)
    {

        for (int i = 0; i < initial_number_of_cells; i++)
        {
            cell cell_;
            
            cell_.x = GenerateRandomNumber(radius + x_min, x_max - radius);
            cell_.y = GenerateRandomNumber(radius + y_min, y_max - radius);
            cell_.vx = 0.0;
            cell_.vy = 0.0;
            cell_.fx.resize(number_of_threads,0.0);
            cell_.fy.resize(number_of_threads,0.0);
            cell_.k = GenerateRandomNumber(k_min, k_max);
            cell_.damp_ratio = GenerateRandomNumber(damping_min, damping_max);
            cell_.mass = 1.0;
            cell_.is_wall = -1;
            
            cells.cells_collection.push_back(cell_);
            
            cells.number_of_cells++;
        }
    }

    void add_wall(
        double x1_, double y1_,
        double x2_, double y2_,
        int number_of_cells_, double youngs_modulus, double damping_ratio, double mass,
        int number_of_threads)
    {
        wall_ids.push_back(number_of_walls);

        double dist_x = x2_ - x1_;
        double dist_y = y2_ - y1_;

        double distance = sqrt(dist_x * dist_x + dist_y * dist_y);

        double dx = dist_x / (double)number_of_cells_;
        double dy = dist_y / (double)number_of_cells_;

        for (int j = 0; j < number_of_cells_; j++)
        {
            cell cell_;
            cell_.x = x1_ + dx * (double)j;
            cell_.y = y1_ + dy * (double)j;
            cell_.vx = 0.0;
            cell_.vy = 0.0;
            cell_.fx.resize(number_of_threads,0.0);
            cell_.fy.resize(number_of_threads,0.0);
            cell_.k=youngs_modulus;
            cell_.damp_ratio = damping_ratio;
            cell_.mass = mass;
            cell_.is_wall = number_of_walls;
            
            cells.cells_collection.push_back(cell_);
            cells.number_of_cells++;
        }
        number_of_walls += 1;
    }

    void start_simulation(double output_interval, double timestep, double duration,int number_of_threads, bool print_output)
    {
        double time_finding_contacts = 0.0;
        double time_in_forces = 0.0;
        double time_in_move = 0.0;
        int n_timesteps=0, n_searchs=0;
        int number_of_cells = cells.number_of_cells;

        omp_set_num_threads(number_of_threads);

        double max_velocity_search = radius * search_multiplier / timestep;

        double simulation_time = 0.0;

        // visualization stuff

        if (print_output)
        {
            double next_output = 0.0;

            const unsigned char bluegreen[] = {0, 170, 255};
            const unsigned char black[] = {0, 0, 0};

            CImg<unsigned char> bg((int)width, (int)height, 1, 3, 255);

            bg.draw_rectangle(0, 0, (int)width, (int)height, bluegreen);

            CImgDisplay dsp((int)width, (int)height, "EvolutionSimulator", 0);

            dsp.display(bg);

            srand((unsigned int)time(NULL)); // randomize seed

            CImg<unsigned char> img(bg);

            Sleep(3000);
            while (!dsp.is_closed() && !dsp.is_keyESC())
            {
                img = bg;

                // display only for output intervals
                if (simulation_time >= next_output)
                {
                    for (int i = 0; i < number_of_cells; i++)
                    {
                        img.draw_circle((int)width - (int)cells.cells_collection[i].x, (int)height - (int)cells.cells_collection[i].y, (int)radius, black);
                    }
                    dsp.display(img);
                    next_output += output_interval;
                }

                do_one_iteration(cells, grid, radius, timestep,
                                 width, height, friction_coefficient, gravity, has_boundaries,
                                 time_finding_contacts,
                                 time_in_forces,
                                 time_in_move, n_timesteps, n_searchs,max_velocity_search,number_of_threads);

                simulation_time += timestep;
                if (simulation_time >= duration)
                    break;
            }
        }
        else
        {
            while (1)
            {

                do_one_iteration(cells, grid, radius,
                                 timestep, width, height, friction_coefficient,
                                 gravity, has_boundaries,
                                 time_finding_contacts,
                                 time_in_forces,
                                 time_in_move, n_timesteps, n_searchs,max_velocity_search,number_of_threads);

                simulation_time += timestep;
                if (simulation_time >= duration)
                    break;
            }
        }
        std::cout << "time_finding_contacts: " << time_finding_contacts << "\n";
        std::cout << "time_in_forces: " << time_in_forces << "\n";
        std::cout << "time_in_move: " << time_in_move << "\n";
        std::cout << "n_timesteps: " << n_timesteps << "\n";
        std::cout << "n_searchs: " << n_searchs << "\n";
    }
};
