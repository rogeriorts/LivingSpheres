#pragma once
#define _USE_MATH_DEFINES

#include "grid.h"
#include "cells.h"

#define M_PI 3.14159265358979323846

struct Simulator
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
    int number_of_threads;
    double timestep = 1e-5;
    double max_velocity_search = 0.0;
    double time_finding_contacts = 0.0;
    double time_in_forces = 0.0;
    double time_in_move = 0.0;
    int n_timesteps = 0;
    int n_searchs = 0;
    int num_contacts_max = 0;

    std::vector<double> wall_ids;
    int number_of_walls;

    double has_boundaries;

    double search_multiplier;

    Simulator(
        double radius_ = 3.0,
        double width_ = 50.0,
        double height_ = 50.0,
        double grid_size_multiplier_ = 2.0,
        double search_multiplier_ = 1e-5,
        double friction_coefficient_ = 0.3,
        int number_of_threads_ = 1,
        bool has_boundaries_ = true,
        double gravity_ = -9.81,
        int num_contacts_max_ = 6
        )
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
        number_of_threads = number_of_threads_;
        num_contacts_max = num_contacts_max_;
        grid.grid_initialize(radius, width, height, grid_size_multiplier_);
    }

    void initialize_living_cells(int initial_number_of_cells,
                                        double spring_coefficient,
                                        double damping_ratio, 
                                        double x_min, double y_min,
                                        double x_max, double y_max, double density,
                                        double temperature,
                                        double heat_conduction_coefficient,
                                        double specific_heat,
                                        bool fixed = false
                                        );

    void add_wall(
        double x1_, double y1_,
        double x2_, double y2_,
        int number_of_cells_, double spring_coefficient,
        double damping_ratio, double density,
                                        double temperature,
                                        double heat_conduction_coefficient,
                                        double specific_heat);

    void start_simulation(double output_interval,
                          double timestep_, double duration, bool print_output);

    bool verify_boundary_contact(double cell_position,
                                 double boundary_max, double &boundary_position);

    bool add_cell_contact();
    
    void add_cells_to_grid();

    void find_contact_between_pair_of_cells(int i, int j);

    void find_contacts_grid();

    void find_contacts_grid_cuda();

    void reset_forces();

    void calculate_force_due_contacts();

    void compute_new_position();

    bool search_trigger();

    void search_step(bool force_search = false);

    void simulation_setup(double timestep_);

    void do_one_iteration();
};