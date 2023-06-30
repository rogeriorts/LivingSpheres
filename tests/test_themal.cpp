#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../include/LivingSpheres/simulator.h"
using namespace std;

TEST_CASE("Test thermal calculation")
{

        // simulation settings
    bool print_animation = false;
    int number_of_threads = 2;
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
    
    //double grid parameters
    double grid_multiplier = 2.0;
    double search_multiplier = 0.0000001;

    // initialize environment
    Simulator sim = Simulator(
        radius, width, height, grid_multiplier,
        search_multiplier, friction_coefficient, number_of_threads,
        true, gravity);

    bool fixed = true;
    double cells_temperature = 300;
    double wall_temperature = 400;

    sim.initialize_living_cells(
        number_of_cells_alive, 1e3, 0.4, 10.0, 200.0, 
        width, height, density,
        cells_temperature,heat_conduction_coefficient,specific_heat,fixed);

    sim.cells.cells_collection[0].x = 25.0;
    sim.cells.cells_collection[0].y = 25.0 + 1.5 * radius;

    sim.add_wall(25.0, 25.0, 25.0, 25.000000001,1, 1e4, 0.1, 
        density,wall_temperature,heat_conduction_coefficient,specific_heat);

    sim.start_simulation(1,timestep,timestep,false);

    CHECK(sim.cells.cells_collection[0].fixed == true);
    CHECK(sim.cells.cells_collection[1].fixed == true);

    CHECK(sim.cells.cells_collection[0].x == Approx(25.0).margin(1e-5));
    CHECK(sim.cells.cells_collection[0].y == Approx(25.0 + 1.5 * radius).margin(1e-5));
    
    CHECK(sim.cells.cells_collection[1].x == Approx(25.0).margin(1e-5));
    CHECK(sim.cells.cells_collection[1].y == Approx(25.0).margin(1e-5));

    CHECK(sim.cells.cells_collection[1].temperature == Approx(400).margin(1e-5));

    double dist_x = sim.cells.cells_collection[1].x - sim.cells.cells_collection[0].x;
    double dist_y = sim.cells.cells_collection[1].y - sim.cells.cells_collection[0].y;

    double distance = sqrt(dist_x*dist_x + dist_y*dist_y);

    double dist_to_contact = distance / 2.0;

    double squared_half_chord = radius*radius - dist_to_contact*dist_to_contact;

    double contact_area = squared_half_chord * M_PI;

    double heat_transfer_rate = heat_conduction_coefficient * (wall_temperature - cells_temperature) * contact_area / distance;

    double calc_heat_transfer_rate = sim.cells.cells_collection[0].transferred_heat[0] + sim.cells.cells_collection[0].transferred_heat[1];

    CHECK(heat_transfer_rate == Approx(calc_heat_transfer_rate).margin(1e-5));

    double new_temperature = cells_temperature + heat_transfer_rate * timestep / 
        (sim.cells.cells_collection[0].mass * sim.cells.cells_collection[0].specific_heat);

    CHECK(new_temperature == Approx(sim.cells.cells_collection[0].temperature).margin(1e-5));
    CHECK(wall_temperature == Approx(sim.cells.cells_collection[1].temperature).margin(1e-5));


}

TEST_CASE("Test thermal vs analytical solution")
{

        // simulation settings
    bool print_animation = false;
    int number_of_threads = 2;
    double timestep = 0.5e-2;
    double duration = 100.0 * timestep;

    // environment settings
    int number_of_cells_alive = 1;
    double gravity = -9.81;
    double width = 50.0;
    double height = 50.0;

    // cell and material settings
    double radius = 3.0; // m
    double friction_coefficient = 0.3;
    double density = 1.0 / (4.0 / 3.0 * M_PI * radius * radius * radius); //kg / m3
    double heat_conduction_coefficient = 0.187*1000; //  W/(m K) ;
    double specific_heat = 1.67e3; //J/(kg K)
    
    //double grid parameters
    double grid_multiplier = 2.0;
    double search_multiplier = 0.00001;

    // initialize environment
    Simulator sim = Simulator(
        radius, width, height, grid_multiplier,
        search_multiplier, friction_coefficient, number_of_threads,
        true, gravity);

    bool fixed = true;
    double cells_temperature = 300;
    double wall_temperature = 400;

    sim.initialize_living_cells(
        number_of_cells_alive, 1e3, 0.4, 10.0, 200.0, 
        width, height, density,
        cells_temperature,heat_conduction_coefficient,specific_heat,fixed);

    sim.cells.cells_collection[0].x = 25.0;
    sim.cells.cells_collection[0].y = 25.0 + 1.5 * radius;

    sim.add_wall(25.0, 25.0, 25.0, 25.000000001,1, 1e4, 0.1, 
        density,wall_temperature,heat_conduction_coefficient,specific_heat);

    sim.start_simulation(1,timestep,duration,false);

    CHECK(sim.cells.cells_collection[0].fixed == true);
    CHECK(sim.cells.cells_collection[1].fixed == true);

    CHECK(sim.cells.cells_collection[0].x == Approx(25.0).margin(1e-5));
    CHECK(sim.cells.cells_collection[0].y == Approx(25.0 + 1.5 * radius).margin(1e-5));
    
    CHECK(sim.cells.cells_collection[1].x == Approx(25.0).margin(1e-5));
    CHECK(sim.cells.cells_collection[1].y == Approx(25.0).margin(1e-5));

    CHECK(sim.cells.cells_collection[1].temperature == Approx(400).margin(1e-5));

    double dist_x = sim.cells.cells_collection[1].x - sim.cells.cells_collection[0].x;
    double dist_y = sim.cells.cells_collection[1].y - sim.cells.cells_collection[0].y;

    double distance = sqrt(dist_x*dist_x + dist_y*dist_y);

    double dist_to_contact = distance / 2.0;

    double squared_half_chord = radius*radius - dist_to_contact*dist_to_contact;

    double contact_area = squared_half_chord * M_PI;

    double final_temperature = 
        wall_temperature - (wall_temperature - cells_temperature) *
        exp( - heat_conduction_coefficient * contact_area * duration / 
        (sim.cells.cells_collection[0].mass * specific_heat * distance)
        );

    CHECK(final_temperature == Approx(sim.cells.cells_collection[0].temperature).margin(0.1));


}
