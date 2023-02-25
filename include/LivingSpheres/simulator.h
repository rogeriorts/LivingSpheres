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

#define M_PI 3.14159265358979323846

using namespace cimg_library;

struct cells_type
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> fx;
    std::vector<double> fy;
    std::vector<double> k;
    std::vector<double> mass;
    std::vector<double> damp_ratio;
    std::vector<int> is_wall;
    int number_of_cells = 0;
};

struct cell_contacts
{
    std::vector<int> cell_a;
    std::vector<int> cell_b;
    std::vector<double> contact_vector_x;
    std::vector<double> contact_vector_y;
    int number_of_contacts = 0;
    int max_number_of_contacts = 0;

    void add_contact(int cell_a_, int cell_b_, double contact_vector_x_,
                     double contact_vector_y_, int last_contact_index)
    {
        if (last_contact_index >= cell_a.size())
        {
            cell_a.push_back(cell_a_);
            cell_b.push_back(cell_b_);
            contact_vector_x.push_back(contact_vector_x_);
            contact_vector_y.push_back(contact_vector_y_);
        }
        else
        {
            cell_a[last_contact_index] = cell_a_;
            cell_b[last_contact_index] = cell_b_;
            contact_vector_x[last_contact_index] = contact_vector_x_;
            contact_vector_y[last_contact_index] = contact_vector_y_;
        }
        number_of_contacts = last_contact_index + 1;
    }
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

void find_contacts(cells_type &cells, cell_contacts &contacts, double radius, bool has_boundaries,
                   double width = 500.0, double height = 500.0)
{

    int contact_count = 0;
    contacts.number_of_contacts = 0;
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        for (int j = i + 1; j < cells.number_of_cells; j++)
        {
            if (cells.is_wall[j] > -1 && cells.is_wall[i] > -1)
                continue;

            double dist_x = cells.x[j] - cells.x[i];
            double dist_y = cells.y[j] - cells.y[i];
            double squared_distance = dist_x * dist_x + dist_y * dist_y;

            if (squared_distance < 4.0 * radius * radius)
            {
                contacts.add_contact(i, j, dist_x, dist_y, contact_count);
                contact_count++;
            }
        }
        if (has_boundaries)
        {
            double boundary_position;
            // verify horizontal contacts with boundaries
            if (verify_boundary_contact(cells.x[i], width, radius, boundary_position))
            {
                double dist_x = boundary_position - cells.x[i];
                contacts.add_contact(i, -1, dist_x, 0.0, contact_count);
                contact_count++;
            }
            if (verify_boundary_contact(cells.y[i], height, radius, boundary_position))
            {
                double dist_y = boundary_position - cells.y[i];
                contacts.add_contact(i, -1, 0.0, dist_y, contact_count);
                contact_count++;
            }
        }
    }
}

void add_cells_to_grid(cells_type &cells, grid_type &grid)
{
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        grid.add_cell_to_block(cells.x[i], cells.y[i], i);
    }
}

void find_contact_between_pair_of_cells(
    cells_type &cells, cell_contacts &contacts,
    int &contact_count, int i, int j, double radius)
{
    if (cells.is_wall[i] > -1 && cells.is_wall[j] > -1)
        return;

    double dist_x = cells.x[j] - cells.x[i];
    double dist_y = cells.y[j] - cells.y[i];
    double squared_distance = dist_x * dist_x + dist_y * dist_y;

    if (squared_distance < 4.0 * radius * radius)
    {
        contacts.add_contact(i, j, dist_x, dist_y, contact_count);
        contact_count++;
    }
}

void find_contacts_grid(cells_type &cells, cell_contacts &contacts, grid_type &grid,
                        double radius, bool has_boundaries, double width = 500.0, double height = 500.0)
{
    int contact_count = 0;
    contacts.number_of_contacts = 0;
    for (int grid_i = 0; grid_i < grid.nx; grid_i++)
    {
        for (int grid_j = 0; grid_j < grid.ny; grid_j++)
        {
            int k_a = grid_i + grid.nx * grid_j;
            for (int c_a = 0; c_a < grid.blocks[k_a].next_position_on_block; c_a++)
            {
                int i = grid.blocks[k_a].cell_ids[c_a];
                // looking for the contacts inside the same grid block
                for (int c_b = c_a + 1; c_b < grid.blocks[k_a].next_position_on_block; c_b++)
                {
                    int j = grid.blocks[k_a].cell_ids[c_b];
                    find_contact_between_pair_of_cells(cells, contacts, contact_count, i, j, radius);
                }
                
                // now, lets compare with the next block i + 1
                if (grid_i + 1 < grid.nx)
                {
                    int k_b = k_a + 1;
                    for (int c_b = 0; c_b < grid.blocks[k_b].next_position_on_block; c_b++)
                    {
                        int j = grid.blocks[k_b].cell_ids[c_b];
                        find_contact_between_pair_of_cells(cells, contacts, contact_count, i, j, radius);
                    }
                }
                // now, lets compare with the next block j + 1
                if (grid_j + 1 < grid.ny)
                {
                    int k_b = k_a + grid.nx;
                    for (int c_b = 0; c_b < grid.blocks[k_b].next_position_on_block; c_b++)
                    {
                        int j = grid.blocks[k_b].cell_ids[c_b];
                        find_contact_between_pair_of_cells(cells, contacts, contact_count, i, j, radius);
                    }
                }
                // now, lets compare with the next block i + 1 and j-1
                if (grid_j - 1 >= 0 && grid_i + 1 < grid.nx)
                {
                    int k_b = k_a - grid.nx + 1;
                    for (int c_b = 0; c_b < grid.blocks[k_b].next_position_on_block; c_b++)
                    {
                        int j = grid.blocks[k_b].cell_ids[c_b];
                        find_contact_between_pair_of_cells(cells, contacts, contact_count, i, j, radius);
                    }
                }
                // now, lets compare with the next block i + 1 and j+1
                if (grid_j + 1 < grid.ny && grid_i + 1 < grid.nx)
                {
                    int k_b = k_a + grid.nx + 1;
                    for (int c_b = 0; c_b < grid.blocks[k_b].next_position_on_block; c_b++)
                    {
                        int j = grid.blocks[k_b].cell_ids[c_b];
                        find_contact_between_pair_of_cells(cells, contacts, contact_count, i, j, radius);
                    }
                }
                
                if (has_boundaries)
                {
                    double boundary_position;
                    // verify horizontal contacts with boundaries
                    if (verify_boundary_contact(cells.x[i], width, radius, boundary_position))
                    {
                        double dist_x = boundary_position - cells.x[i];
                        contacts.add_contact(i, -1, dist_x, 0.0, contact_count);
                        contact_count++;
                    }
                    if (verify_boundary_contact(cells.y[i], height, radius, boundary_position))
                    {
                        double dist_y = boundary_position - cells.y[i];
                        contacts.add_contact(i, -1, 0.0, dist_y, contact_count);
                        contact_count++;
                    }
                }
            }
        }
    }
}

void reset_forces(cells_type &cells)
{
    fill(cells.fx.begin(), cells.fx.end(), 0.0);
    fill(cells.fy.begin(), cells.fy.end(), 0.0);
}

void calculate_force_due_contacts(cells_type &cells, cell_contacts &contacts, double radius, double friction_coefficient)
{

    for (int c = 0; c < contacts.number_of_contacts; c++)
    {
        double dist_x = contacts.contact_vector_x[c];
        double dist_y = contacts.contact_vector_y[c];
        double squared_distance = dist_x * dist_x + dist_y * dist_y;

        // calculate normal vector between the contacting cells
        double distance = sqrt(squared_distance);
        double x_normal = dist_x / distance;
        double y_normal = dist_y / distance;

        double x_tangential = y_normal;
        double y_tangential = -x_normal;

        int i = contacts.cell_a[c];
        int j = contacts.cell_b[c];

        double avg_k, avg_mass, avg_damp_ratio, vx_b, vy_b, spring_deformation;

        if (j < 0)
        {
            avg_k = cells.k[i];
            avg_mass = cells.mass[i];
            avg_damp_ratio = cells.damp_ratio[i];
            vx_b = 0.0;
            vy_b = 0.0;
            spring_deformation = radius - distance;
        }
        else
        {
            avg_k = (cells.k[i] + cells.k[j]) * 0.5;
            avg_mass = (cells.mass[i] + cells.mass[j]) * 0.5;
            avg_damp_ratio = (cells.damp_ratio[i] + cells.damp_ratio[j]) * 0.5;
            vx_b = cells.vx[j];
            vy_b = cells.vy[j];
            spring_deformation = 2.0 * radius - distance;
        }
        // spring force
        double spring_force = std::min(avg_k * spring_deformation, avg_k * radius);
        // damping force

        double vx_rel = vx_b - cells.vx[i];
        double vy_rel = vy_b - cells.vy[i];

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

        cells.fx[i] += fx;
        cells.fy[i] += fy;
        if (j >= 0)
        {
            cells.fx[j] -= fx;
            cells.fy[j] -= fy;
        }
    }
}

void compute_new_position(cells_type &cells, double timestep, double gravity)
{
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        // gravity influence
        if (cells.is_wall[i] < 0)
        {
            // calculating accelerations
            double ax = cells.fx[i] / cells.mass[i];
            double ay = cells.fy[i] / cells.mass[i] + gravity;

            // calculate new positions
            //  s = s0 + v0 * t + a * t^2 /2
            double acceleration_timeste_factor = timestep * timestep * 0.5;
            cells.x[i] += cells.vx[i] * timestep + ax * acceleration_timeste_factor;
            cells.y[i] += cells.vy[i] * timestep + ay * acceleration_timeste_factor;

            // calculating new velocities
            // v = v0 + a * t
            cells.vx[i] += ax * timestep;
            cells.vy[i] += ay * timestep;
        }
    }
}

using namespace std::chrono;

void do_one_iteration(cells_type &cells, cell_contacts &contacts,grid_type &grid,
                      double radius, double timestep, double width, double height, double friction_coefficient, double gravity,
                      bool has_boundaries, bool use_grid,
                      double &time_finding_contacts,
                    double &time_in_forces, double &time_in_move)
{
    auto start = high_resolution_clock::now();
    if(use_grid)
    {
        
        grid.reset_grid();
        add_cells_to_grid(cells, grid);
        find_contacts_grid(cells,contacts,grid,radius,has_boundaries,width,height);
    }
    else{
        find_contacts(cells, contacts, radius, has_boundaries, width, height);
    }
    auto stop = high_resolution_clock::now();
    time_finding_contacts += duration<double>(stop - start).count();
    
    
    start = high_resolution_clock::now();
    reset_forces(cells);
    calculate_force_due_contacts(cells, contacts, radius, friction_coefficient);
    stop = high_resolution_clock::now();
    time_in_forces += duration<double>(stop - start).count();
    // if(!has_boundaries) calculate_force_due_boundaries(radius, cells, timestep, width, height, friction_coefficient);
    start = high_resolution_clock::now();
    compute_new_position(cells, timestep, gravity);
    stop = high_resolution_clock::now();
    time_in_move += duration<double>(stop - start).count();
}

struct Environment
{
    cells_type cells;
    cell_contacts contacts;
    grid_type grid;

    double radius;
    double friction_coefficient;
    double initial_number_of_cells;
    double width;
    double height;
    double gravity;
    double minimum_timesteps_per_loading;
    bool cells_initialized;
    bool walls_initialized;

    std::vector<double> wall_ids;
    int number_of_walls;

    double has_boundaries;

    double timestep;

    Environment(
        double radius_,
        double width_,
        double height_,
        double grid_multiplier,
        double friction_coefficient_,
        bool has_boundaries_ = true,
        double gravity_ = -9.81,
        double minimum_timesteps_per_loading_ = 15)
    {
        radius = radius_;
        height = height_;
        width = width_;
        friction_coefficient = friction_coefficient_;
        gravity = gravity_;
        has_boundaries = has_boundaries_;
        minimum_timesteps_per_loading = minimum_timesteps_per_loading_;
        cells_initialized = false;
        walls_initialized = false;
        number_of_walls = 0;
        grid.grid_initialize(radius, width, height,grid_multiplier);
    }

    void initialize_living_cells(int initial_number_of_cells,
                                 double k_min, double k_max,
                                 double damping_min, double damping_max, double x_min, double y_min,
                                 double x_max, double y_max)
    {

        for (int i = 0; i < initial_number_of_cells; i++)
        {
            cells.x.push_back(GenerateRandomNumber(radius + x_min, x_max - radius));
            cells.y.push_back(GenerateRandomNumber(radius + y_min, y_max - radius));
            cells.vx.push_back(0.0);
            cells.vy.push_back(0.0);
            cells.fx.push_back(0.0);
            cells.fy.push_back(0.0);
            cells.k.push_back(GenerateRandomNumber(k_min, k_max));
            cells.damp_ratio.push_back(GenerateRandomNumber(damping_min, damping_max));
            cells.mass.push_back(1.0);
            cells.is_wall.push_back(-1);
            cells.number_of_cells++;
        }
    }

    void add_wall(
        double x1_, double y1_,
        double x2_, double y2_,
        int number_of_cells_, double youngs_modulus, double damping_ratio, double mass)
    {
        wall_ids.push_back(number_of_walls);

        double dist_x = x2_ - x1_;
        double dist_y = y2_ - y1_;

        double distance = sqrt(dist_x * dist_x + dist_y * dist_y);

        double dx = dist_x / (double)number_of_cells_;
        double dy = dist_y / (double)number_of_cells_;

        for (int j = 0; j < number_of_cells_; j++)
        {
            cells.x.push_back(x1_ + dx * (double)j);
            cells.y.push_back(y1_ + dy * (double)j);
            cells.vx.push_back(0.0);
            cells.vy.push_back(0.0);
            cells.fx.push_back(0.0);
            cells.fy.push_back(0.0);
            cells.k.push_back(youngs_modulus);
            cells.damp_ratio.push_back(damping_ratio);
            cells.mass.push_back(mass);
            cells.is_wall.push_back(number_of_walls);
            cells.number_of_cells++;
        }
        number_of_walls += 1;
    }

    void start_simulation(double output_interval, double duration,bool use_grid, bool print_output)
    {
        double time_finding_contacts = 0.0;
        double time_in_forces = 0.0;
        double time_in_move = 0.0;
        int number_of_cells = cells.number_of_cells;

        // CALCULATE TIMESTEP
        double max_k = cells.k[0];
        double min_mass = cells.mass[0];

        // this criteria is not the best one, someday I'll improve it
        for (int i = 0; i < number_of_cells; i++)
        {
            if (cells.k[i] > max_k)
                max_k = cells.k[i];
            if (cells.mass[i] < min_mass)
                min_mass = cells.mass[i];
        }

        double cell_stiffness = max_k * 2.0 * (double)radius;

        double timestep = M_PI / (2 * (double)minimum_timesteps_per_loading) * sqrt(min_mass / cell_stiffness);

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
                        img.draw_circle((int)width - (int)cells.x[i], (int)height - (int)cells.y[i], (int)radius, black);
                    }
                    dsp.display(img);
                    next_output += output_interval;
                    Sleep(100);
                }

                do_one_iteration(cells, contacts,grid, radius, timestep, 
                width, height, friction_coefficient, gravity, has_boundaries, use_grid,
                    time_finding_contacts,
                    time_in_forces,
                    time_in_move
                );

                simulation_time += timestep;
                if (simulation_time >= duration)
                    break;
            }
        }
        else
        {
            while (1)
            {

                do_one_iteration(cells, contacts,grid, radius, 
                timestep, width, height, friction_coefficient, 
                gravity, has_boundaries, use_grid,
                    time_finding_contacts,
                    time_in_forces,
                    time_in_move
                );

                simulation_time += timestep;
                if (simulation_time >= duration)
                    break;
            }
        }
        std::cout<< "time_finding_contacts" <<time_finding_contacts<<"\n";
        std::cout<< "time_in_forces" <<time_in_forces<<"\n";
        std::cout<< "time_in_move" <<time_in_move<<"\n";
    }
};
