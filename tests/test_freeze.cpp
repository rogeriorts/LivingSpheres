#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../include/LivingSpheres/simulator.h"
using namespace std;

TEST_CASE("Test freeze cells")
{
    // simulation settings
    bool print_animation = false;
    int number_of_threads = 1;
    double timestep = 0.5e-2;

    // environment settings
    int number_of_cells_alive = 1;
    double gravity = -9.81;
    double width = 50.0;
    double height = 50.0;

    // cell and material settings
    double radius = 3.0; // m
    double friction_coefficient = 0.3;
    double density = 1.0 / (4.0 / 3.0 * M_PI * radius * radius * radius); //kg / m3
    double heat_conduction_coefficient = 0.187; //  W/(m K) ;
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

    bool fixed = true;

    sim.initialize_living_cells(
        number_of_cells_alive, 1e3, 0.4, 10.0, 200.0, 
        width, height, density,
        cells_temperature,heat_conduction_coefficient,specific_heat,fixed);

    sim.cells.cells_collection[0].x = width / 2.0 + 4*radius;
    sim.cells.cells_collection[0].y = height / 2.0;

    sim.initialize_living_cells(
        number_of_cells_alive, 1e3, 0.4, 10.0, 200.0, 
        width, height, density,
        cells_temperature,heat_conduction_coefficient,specific_heat);

    sim.cells.cells_collection[1].x = width / 2.0 - 4*radius;
    sim.cells.cells_collection[1].y = height / 2.0;

    CHECK(sim.cells.cells_collection[0].fixed == true);
    CHECK(sim.cells.cells_collection[1].fixed == false);

    sim.do_one_iteration();

    CHECK(sim.cells.cells_collection[0].x == Approx(width / 2.0 + 4*radius).margin(1e-5));
    CHECK(sim.cells.cells_collection[0].y == Approx(height / 2.0).margin(1e-5));
    
    CHECK(sim.cells.cells_collection[1].x == Approx(width / 2.0 - 4*radius).margin(1e-5));
    CHECK(sim.cells.cells_collection[1].y == Approx(
        height / 2.0 + gravity * timestep * timestep * 0.5 ).margin(1e-5));

}

