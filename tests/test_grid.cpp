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
    cells_type cells;
    cell_contacts contacts;

    double width = 40.0;
    double height = 100.0;

    double radius = 4.24;
    grid_type grid;
    grid.grid_initialize(radius,width,height);
    double gravity = -9.81;

    cells.x.push_back(25.0);
    cells.y.push_back(95.0);
    cells.vx.push_back(3.0);
    cells.vy.push_back(0.0);
    cells.fx.push_back(0.0);
    cells.fy.push_back(0.0);
    cells.k.push_back(1.0e6);
    cells.damp_ratio.push_back(0.1);
    cells.mass.push_back(1.0);
    cells.is_wall.push_back(-1);
    cells.number_of_cells++;

    grid.add_cell_to_block(cells.x[0],cells.y[0],0);

    CHECK(grid.blocks[19].next_position_on_block == 1);
    CHECK(grid.blocks[19].cell_ids[0] == 0);

    double timestep = 3.0;

    compute_new_position(cells, timestep, gravity);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * 3.0 = 26.5
    CHECK(cells.x[0] == Approx(34.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * 3.0 + (-9.81) * 3.0 * 3.0 * 0.5 = -4.905e-08
    CHECK(cells.y[0] == Approx(50.855000000000004).margin(1e-10));
    // same velocity
    CHECK(cells.vx[0] == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 3
    CHECK(cells.vy[0] == Approx(-29.43).margin(1e-8));

    grid.reset_grid();

    grid.add_cell_to_block(cells.x[0],cells.y[0],0);

    CHECK(grid.blocks[19].next_position_on_block == 0);

    CHECK(grid.blocks[11].next_position_on_block == 1);
    CHECK(grid.blocks[11].cell_ids[0] == 0);

    compute_new_position(cells, timestep, gravity);

    grid.reset_grid();

    CHECK(grid.add_cell_to_block(cells.x[0],cells.y[0],0) ==false);
    for(int i=0 ; i< grid.blocks.size();i++)
    {
        CHECK(grid.blocks[i].next_position_on_block == 0);
    }

}