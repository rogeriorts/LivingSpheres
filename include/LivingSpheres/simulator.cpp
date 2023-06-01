
#include <iostream>
#include "CImg.h"
#include "util.h"
#include <algorithm>
#include <chrono>
#include <omp.h>
#include "simulator.h"

#define M_PI 3.14159265358979323846

bool Simulator::verify_boundary_contact(double cell_position, double boundary_max,
                                        double &boundary_position)
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

void Simulator::add_cells_to_grid()
{
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        cells.cells_collection[i].grid_block = grid.add_cell_to_block(
            cells.cells_collection[i].x, cells.cells_collection[i].y, i);
    }
}

void Simulator::find_contact_between_pair_of_cells(int i, int j)
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

void Simulator::find_contacts_grid()
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

        for (int grid_i_neighbor = std::max(0, grid_i - 1);
             grid_i_neighbor <= std::min(grid.nx - 1, grid_i + 1);
             grid_i_neighbor++)
        {
            for (int grid_j_neighbor = std::max(0, grid_j - 1);
                 grid_j_neighbor <= std::min(grid.ny - 1, grid_j + 1);
                 grid_j_neighbor++)
            {
                int k_neighbor = grid_i_neighbor + grid.nx * grid_j_neighbor;
                for (int c_b = 0; c_b < grid.blocks[k_neighbor].next_position_on_block; c_b++)
                {
                    int j = grid.blocks[k_neighbor].cell_ids[c_b];
                    find_contact_between_pair_of_cells(i, j);
                }
            }
        }

        if (has_boundaries)
        {
            double boundary_position;
            // verify horizontal contacts with boundaries
            if (verify_boundary_contact(cells.cells_collection[i].x, width, boundary_position))
            {
                double dist_x = boundary_position - cells.cells_collection[i].x;
                cells.cells_collection[i].add_contact(-1, dist_x, 0.0);
            }
            if (verify_boundary_contact(cells.cells_collection[i].y, height, boundary_position))
            {
                double dist_y = boundary_position - cells.cells_collection[i].y;
                cells.cells_collection[i].add_contact(-1, 0.0, dist_y);
            }
        }
    }
}

void Simulator::reset_forces()
{
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        std::fill(cells.cells_collection[i].fx.begin(),
                  cells.cells_collection[i].fx.end(), 0.0);

        std::fill(cells.cells_collection[i].fy.begin(),
                  cells.cells_collection[i].fy.end(), 0.0);
    }
}

void Simulator::calculate_force_due_contacts()
{
#pragma omp parallel for
    for (int cell = 0; cell < cells.number_of_cells; cell++)
    {
        int t = omp_get_thread_num();

        for (int contact = 0; contact < cells.cells_collection[cell].number_of_contacts; contact++)
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
                avg_k = cells.cells_collection[i].spring_coefficient;
                avg_mass = cells.cells_collection[i].mass;
                avg_damp_ratio = cells.cells_collection[i].damping_ratio;
                vx_b = 0.0;
                vy_b = 0.0;
                spring_deformation = radius - distance;
            }
            else
            {
                avg_k = (cells.cells_collection[i].spring_coefficient + cells.cells_collection[j].spring_coefficient) * 0.5;
                avg_mass = (cells.cells_collection[i].mass + cells.cells_collection[j].mass) * 0.5;
                avg_damp_ratio = (cells.cells_collection[i].damping_ratio + cells.cells_collection[j].damping_ratio) * 0.5;
                vx_b = cells.cells_collection[j].vx;
                vy_b = cells.cells_collection[j].vy;
                spring_deformation = 2.0 * radius - distance;
            }
            // spring force
            if (spring_deformation < 0.0)
                continue;

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

void Simulator::compute_new_position()
{
#pragma omp parallel for
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        // gravity influence
        if (cells.cells_collection[i].is_wall < 0)
        {
            // calculating accelerations
            double fx = 0.0;
            double fy = 0.0;

            for (int t = 0; t < number_of_threads; t++)
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

bool Simulator::search_trigger()
{

    double squared_max_velocity = max_velocity_search * max_velocity_search;
    for (int i = 0; i < cells.number_of_cells; i++)
    {
        if (cells.cells_collection[i].vx * cells.cells_collection[i].vx +
                cells.cells_collection[i].vy * cells.cells_collection[i].vy >
            squared_max_velocity)
            return true;
    }
    return false;
}

void Simulator::do_one_iteration()
{
    auto start = std::chrono::high_resolution_clock::now();

    if (search_trigger())
    {

        grid.reset_grid();
        add_cells_to_grid();
        find_contacts_grid();
        n_searchs++;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    time_finding_contacts += std::chrono::duration<double>(stop - start).count();

    start = std::chrono::high_resolution_clock::now();

    reset_forces();
    calculate_force_due_contacts();

    stop = std::chrono::high_resolution_clock::now();
    time_in_forces += std::chrono::duration<double>(stop - start).count();

    start = std::chrono::high_resolution_clock::now();

    compute_new_position();

    stop = std::chrono::high_resolution_clock::now();
    time_in_move += std::chrono::duration<double>(stop - start).count();
    n_timesteps++;
}

void Simulator::initialize_living_cells(int initial_number_of_cells,
                                        double k_min, double k_max,
                                        double damping_min, double damping_max, double x_min, double y_min,
                                        double x_max, double y_max)
{

    for (int i = 0; i < initial_number_of_cells; i++)
    {
        cell cell_;

        cell_.x = GenerateRandomNumber(radius + x_min, x_max - radius);
        cell_.y = GenerateRandomNumber(radius + y_min, y_max - radius);
        cell_.vx = 0.0;
        cell_.vy = 0.0;
        cell_.fx.resize(number_of_threads, 0.0);
        cell_.fy.resize(number_of_threads, 0.0);
        cell_.spring_coefficient = GenerateRandomNumber(k_min, k_max);
        cell_.damping_ratio = GenerateRandomNumber(damping_min, damping_max);
        cell_.mass = 1.0;
        cell_.is_wall = -1;

        cells.cells_collection.push_back(cell_);

        cells.number_of_cells++;
    }
}

void Simulator::add_wall(
    double x1_, double y1_,
    double x2_, double y2_,
    int number_of_cells_, double spring_coefficient, double damping_ratio, double mass)
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
        cell_.fx.resize(number_of_threads, 0.0);
        cell_.fy.resize(number_of_threads, 0.0);
        cell_.spring_coefficient = spring_coefficient;
        cell_.damping_ratio = damping_ratio;
        cell_.mass = mass;
        cell_.is_wall = number_of_walls;

        cells.cells_collection.push_back(cell_);
        cells.number_of_cells++;
    }
    number_of_walls += 1;
}

void Simulator::start_simulation(double output_interval, double timestep_, double duration, bool print_output)
{

    int number_of_cells = cells.number_of_cells;

    omp_set_num_threads(number_of_threads);

    timestep = timestep_;

    max_velocity_search = radius * search_multiplier / timestep;

    double simulation_time = 0.0;

    srand((unsigned int)time(NULL)); // randomize seed

    // visualization stuff
    if (!print_output)
    {
        while (1)
        {

            do_one_iteration();

            simulation_time += timestep;
            if (simulation_time >= duration)
                break;
        }
    }
    else
    {
        double next_output = 0.0;

        const unsigned char bluegreen[] = {0, 170, 255};
        const unsigned char black[] = {0, 0, 0};

        cimg_library::CImg<unsigned char> bg((int)width, (int)height, 1, 3, 255);

        bg.draw_rectangle(0, 0, (int)width, (int)height, bluegreen);

        cimg_library::CImgDisplay dsp((int)width, (int)height, "EvolutionSimulator", 0);

        dsp.display(bg);

        cimg_library::CImg<unsigned char> img(bg);

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

            do_one_iteration();

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