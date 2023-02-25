#define _USE_MATH_DEFINES
#include <iostream>
#include "../include/LivingSpheres/simulator.h"


int main()
{

    double radius = 3;

    int number_of_cells_alive = 200;
    double output_interval = 1;

    double width = 500.0;
    double height = 500.0;
    double friction_coefficient = 0.3;

    Environment env = Environment(radius, width, height, friction_coefficient,true);

    env.initialize_living_cells(number_of_cells_alive, 1e4, 1.1e4, 0.1, 0.1, 10.0, 200.0, width, height);

    env.add_wall(0.0, 0.0, 500.0, 500.0, 50, 1e4, 0.1,1.0);

    env.add_wall(500.0,0.0,0.0 , 500.0, 50, 1e4, 0.1,1.0);
    
    env.start_simulation(0.1, 100, true);

    return 0;
}