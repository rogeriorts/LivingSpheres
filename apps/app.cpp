#define _USE_MATH_DEFINES
#include <iostream>
#include "../include/LivingSpheres/simulator.h"


int main()
{

    double radius = 3;

    int number_of_cells_alive = 100;
    double output_interval = 0.1;

    double width = 500.0;
    double height = 500.0;
    double friction_coefficient = 0.3;
    double grid_multiplier = 10.0;

    Environment env = Environment(radius, width, height,grid_multiplier, friction_coefficient,true);

    env.initialize_living_cells(number_of_cells_alive, 1e4, 1.1e4, 0.4, 0.5, 10.0, 200.0, width, height);

    env.add_wall(0.0, 0.0, 500.0, 500.0, 50, 1e4, 0.1,1.0);

    env.add_wall(500.0,0.0,0.0 , 500.0, 50, 1e4, 0.1,1.0);
    bool use_grid = true;
    
    auto start = std::chrono::steady_clock::now();
    env.start_simulation(output_interval, 10, use_grid,true);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Duration [seconds]: " << diff.count() << std::endl;

    return 0;
}