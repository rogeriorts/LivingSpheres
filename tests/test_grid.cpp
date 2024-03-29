#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../include/LivingSpheres/simulator.h"
using namespace std;

TEST_CASE("Test one cell being added to the grid")
{

    double width = 40.0;
    double height = 30.0;

    double radius = 2.45;
    grid_type grid;
    grid.grid_initialize(radius,width,height);

    //testing grid_size

    CHECK(grid.nx == 5);
    CHECK(grid.ny == 4);
    CHECK(grid.block_size_x == Approx(8.0).margin(1e-8));
    CHECK(grid.block_size_y == Approx(7.5).margin(1e-8));
    //testing a cell inside the grid   

    double cell_x = 25.0;
    double cell_y = 17.0;
    grid.add_cell_to_block(cell_x,cell_y,15);

    for(int i=0 ; i< grid.blocks.size();i++)
    {
        if (i == 13){
            CHECK(grid.blocks[i].next_position_on_block == 1);
            CHECK(grid.blocks[i].cell_ids[0] == 15);
        }
        else
        {
            CHECK(grid.blocks[i].cell_ids.size() == 0);
            CHECK(grid.blocks[i].next_position_on_block == 0);
        }
    }
    grid.reset_grid();
    for(int i=0 ; i< grid.blocks.size();i++)
    {
        CHECK(grid.blocks[i].next_position_on_block == 0);
        if (i == 13){
            CHECK(grid.blocks[i].cell_ids.size() == 1);
        }
        else
        {
            CHECK(grid.blocks[i].cell_ids.size() == 0);
        }
    }
}

TEST_CASE("One cell falling due gravity inside the grid")
{

    // Create a 200 cells randomly positioned  

    double width = 40.0;
    double height = 100.0;

    double radius = 4.24;
    double gravity = -9.81;

    Simulator sim = Simulator(
        radius, width, height, 3.0, 
        0.001, 0.3, 1, true, gravity);
  

    cell cell_;
    cell_.x = 25.0;
    cell_.y = 95.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.transferred_heat.push_back(0.0);
    cell_.spring_coefficient = 1.0e6;
    cell_.damping_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.is_wall = -1;
    sim.cells.cells_collection.push_back(cell_);
    sim.cells.number_of_cells++;

    int grid_block = sim.grid.add_cell_to_block(cell_.x,cell_.y,0);

    CHECK(grid_block > -1);

    CHECK(sim.grid.blocks[19].next_position_on_block == 1);
    CHECK(sim.grid.blocks[19].cell_ids[0] == 0);

    
    sim.timestep = 3.0;
    sim.compute_new_position();

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * 3.0 = 26.5
    CHECK(sim.cells.cells_collection[0].x == Approx(34.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * 3.0 + (-9.81) * 3.0 * 3.0 * 0.5 = -4.905e-08
    CHECK(sim.cells.cells_collection[0].y == Approx(50.855000000000004).margin(1e-10));
    // same velocity
    CHECK(sim.cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 3
    CHECK(sim.cells.cells_collection[0].vy == Approx(-29.43).margin(1e-8));

    sim.grid.reset_grid();

    grid_block = sim.grid.add_cell_to_block(sim.cells.cells_collection[0].x,sim.cells.cells_collection[0].y,0);

    CHECK(grid_block > -1);

    CHECK(sim.grid.blocks[19].next_position_on_block == 0);

    CHECK(sim.grid.blocks[11].next_position_on_block == 1);
    CHECK(sim.grid.blocks[11].cell_ids[0] == 0);

    sim.compute_new_position();

    sim.grid.reset_grid();

    CHECK(sim.grid.add_cell_to_block(sim.cells.cells_collection[0].x,sim.cells.cells_collection[0].y,0) < 0 );
    for(int i=0 ; i< sim.grid.blocks.size();i++)
    {
        CHECK(sim.grid.blocks[i].next_position_on_block == 0);
    }

}