#define _USE_MATH_DEFINES
#include <iostream>
#include "../include/LivingSpheres/simulator.h"
#include <chrono>

int main()
{

    // simulation settings
    bool print_animation = true;
    int number_of_threads = 6;
    double timestep = 0.5e-2;
    double simulation_duration = 300;
    double output_interval = 0.01;

    // environment settings
    int number_of_cells_alive = 1000;
    double gravity = -9.81;
    double width = 800.0;
    double height = 800.0;

    // cell and material settings
    double radius = 3.0; // m
    double friction_coefficient = 0.3;
    double density = 1.0 / (4.0 / 3.0 * M_PI * radius * radius * radius); //kg / m3
    double heat_conduction_coefficient = 0.187*10000; //  W/(m K) ;
    double specific_heat = 1.67e3; //J/(kg K)
    double cells_temperature = 300;
    double wall_temperature = 400;
    
    //double grid parameters
    double grid_multiplier = 2.0;
    double search_multiplier = 0.0000001;

    // initialize environment
    Simulator sim = Simulator(
        radius, width, height, grid_multiplier,
        search_multiplier, friction_coefficient, number_of_threads,
        true, gravity);

    sim.initialize_living_cells(
        number_of_cells_alive, 1e3, 0.4, 4.0*radius, 4.0*radius, 
        width- 4.0*radius, height - 4.0*radius, density,
        cells_temperature,heat_conduction_coefficient,specific_heat);

    // adding wall cells

    sim.add_wall(radius, radius, radius, height - radius, 100, 1e4, 0.1, 
        density, wall_temperature,heat_conduction_coefficient,specific_heat);
    sim.add_wall(width - radius, radius, width - radius, height - radius, 100, 1e4, 0.1, 
        density,wall_temperature,heat_conduction_coefficient,specific_heat);
    sim.add_wall(radius, height - radius, width - radius, height - radius, 100, 1e4, 0.1, 
        density,wall_temperature,heat_conduction_coefficient,specific_heat);
    sim.add_wall(radius, radius, width - radius, radius, 100, 1e4, 0.1, 
        density,wall_temperature,heat_conduction_coefficient,specific_heat);

    auto start = std::chrono::steady_clock::now();

    sim.start_simulation(output_interval, timestep,
                         simulation_duration, print_animation);

    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end - start;

    std::cout << "Duration [seconds]: " << diff.count() << std::endl;

    return 0;
}