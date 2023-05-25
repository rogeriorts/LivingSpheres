#define _USE_MATH_DEFINES
#include <iostream>
#include "../include/LivingSpheres/simulator.h"

int main()
{
    
    //simulation settings
    bool print_animation = true;
    int number_of_threads = 10;
    
    double timestep = 0.5e-2;
    double simulation_duration = 300;
    double output_interval = 0.1;

    //environment settings
    int number_of_cells_alive = 6000;
    double cell_radius = 3.0;
    double gravity = -9.81;
    double width = 1700.0;
    double height = 800.0;
    
    //material settings
    double friction_coefficient = 0.3;
    double grid_multiplier = 2.0;
    double search_multiplier = 0.001;

    //initialize environment
    Environment env = Environment(
        cell_radius, width, height, grid_multiplier, 
        search_multiplier, friction_coefficient, 
        true, gravity);

    env.initialize_living_cells(
        number_of_cells_alive, 1e3, 1.1e3, 0.4, 0.5, 10.0, 200.0, width, height,number_of_threads
        );
    
    //adding wall cells
    env.add_wall(0.0, 0.0, width, height, 500, 1e4, 0.1, 1.0,number_of_threads);

    env.add_wall(width, 0.0, 0.0, height, 500, 1e4, 0.1, 1.0,number_of_threads);

    // env.add_wall(radius, radius, radius, 500.0 - radius, 100, 1e4, 0.1,1.0);
    // env.add_wall(500.0 - radius, radius, 500.0 - radius, 500.0 - radius, 100, 1e4, 0.1,1.0);
    // env.add_wall(radius, 500.0 - radius, 500.0 - radius, 500.0 - radius, 100, 1e4, 0.1,1.0);
    // env.add_wall(radius, radius, 500.0 - radius, radius, 100, 1e4, 0.1,1.0);

    auto start = std::chrono::steady_clock::now();

    env.start_simulation(output_interval,timestep, 
        simulation_duration, number_of_threads, 
        print_animation);

    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end - start;

    std::cout << "Duration [seconds]: " << diff.count() << std::endl;

    return 0;
}