#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../include/LivingSpheres/simulator.h"
using namespace std;

TEST_CASE("Two timesteps of 2 living cells")
{

    // Create a 200 cells randomly positioned
    cells_type cells;
    grid_type grid;
    
    
    double radius = 4.24;

    double gravity = -9.81;

    cell cell_;
    double width = 100.0;
    double height = 100.0;

    grid.grid_initialize(radius,width,height);
    cell_.x = 25.0;
    cell_.y = 25.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.k = 1.0e6;
    cell_.damp_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.is_wall = -1;
    cells.cells_collection.push_back(cell_);
    cells.number_of_cells++;

    cell_.x = 31.0;
    cell_.y = 31.0;
    cell_.vx = 0.0;
    cell_.vy = -3.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.k = 2.0e6;
    cell_.damp_ratio = 0.2;
    cell_.mass = 2.0;
    cell_.is_wall = -1;
    cells.cells_collection.push_back(cell_);
    cells.number_of_cells++;

    add_cells_to_grid(cells, grid);

    double timestep = .3;
    double friction_coefficient = 0.3;

    find_contacts_grid(cells,grid, radius, false,width,height);
    
    reset_forces(cells);
    calculate_force_due_contacts(cells, radius, friction_coefficient);
    // forces must be zero at this point
    CHECK(cells.cells_collection[0].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.cells_collection[0].fy[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.cells_collection[1].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.cells_collection[1].fy[0] == Approx(0.0).margin(1e-10));
    compute_new_position(cells, timestep, gravity,1);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * .3 = 26.5
    CHECK(cells.cells_collection[0].x == Approx(25.9).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * .3 + (-9.81) * .3 * .3 * 0.5
    CHECK(cells.cells_collection[0].y == Approx(24.55855).margin(1e-10));
    // same velocity
    CHECK(cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 1e-4
    CHECK(cells.cells_collection[0].vy == Approx(-2.943).margin(1e-8));

    // same, since there is no vx at this point
    CHECK(cells.cells_collection[1].x == Approx(31.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y =  31.0 + (-3) * .3 + (-9.81) * .3 * .3 * 0.5
    CHECK(cells.cells_collection[1].y == Approx(29.65855).margin(1e-8));
    CHECK(cells.cells_collection[1].vx == Approx(0.0).margin(1e-8));
    // same velocity since there is no acceleration
    // vy = vy_0 + g * t
    // vy = -3 + -9.81 * .3
    CHECK(cells.cells_collection[1].vy == Approx(-5.943).margin(1e-8));

    // NEW TIMESTEP
    find_contacts_grid(cells,grid, radius, false,width,height);
    // WE MUST HAVE ONE CONTACT HERE
    REQUIRE(cells.cells_collection[0].number_of_contacts == 1);
    reset_forces(cells);
    calculate_force_due_contacts(cells, radius, friction_coefficient);

    CHECK(cells.cells_collection[0].fx[0] == Approx(-1347178.0769591301).margin(1e-7));
    CHECK(cells.cells_collection[0].fy[0] == Approx(-1347178.0769591301).margin(1e-7));
    CHECK(cells.cells_collection[1].fx[0] == Approx(1347178.0769591301).margin(1e-7));
    CHECK(cells.cells_collection[1].fy[0] == Approx(1347178.0769591301).margin(1e-7));

    compute_new_position(cells, timestep, gravity,1);

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 25.9 + 3.0 * .3 + (-1347178.0769591301)/1.0 *.3*.3*0.5= 26.5
    CHECK(cells.cells_collection[0].x == Approx(-60596.21346316085).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 24.55855 + -2.943 * .3 + (-9.81 -1347178.0769591301) * .3 * .3 * 0.5
    CHECK(cells.cells_collection[0].y == Approx(-60599.77926316086).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 3.0 + (-1347178.0769591301) * .3
    CHECK(cells.cells_collection[0].vx == Approx(-404150.42308773904).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -2.943 + (-9.81 -1347178.0769591301) * .3
    CHECK(cells.cells_collection[0].vy == Approx(-404159.3090877391).margin(1e-8));

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 31.0 + 0.0 * .3 + (1347178.0769591301)/2.0 *.3*.3*0.5
    CHECK(cells.cells_collection[1].x == Approx(30342.506731580426).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 29.65855 + -5.943 * .3 + (-9.81 +1347178.0769591301/2.0) * .3 * .3 * 0.5
    CHECK(cells.cells_collection[1].y == Approx(30338.940931580422).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 0.0 + (1347178.0769591301/2.0) * .3
    CHECK(cells.cells_collection[1].vx == Approx(202076.71154386952).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -5.943 + (-9.81 +1347178.0769591301)/2 * .3
    CHECK(cells.cells_collection[1].vy == Approx(202069.2970438695).margin(1e-8));
}

TEST_CASE("One cell falling due gravity")
{

    // Create a 200 cells randomly positioned
    cells_type cells;
    grid_type grid;

    double radius = 4.24;

    double gravity = -9.81;

    cell cell_;
    double width = 100.0;
    double height = 100.0;

    grid.grid_initialize(radius,width,height);

    cell_.x = 25.0;
    cell_.y = 25.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.k = 1.0e6;
    cell_.damp_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.is_wall = -1;
    cells.cells_collection.push_back(cell_);
    cells.number_of_cells++;
    
    add_cells_to_grid(cells, grid);

    double timestep = 3.0;
    double friction_coefficient = 0.3;

    find_contacts_grid(cells,grid, radius, false,width,height);
    REQUIRE(cells.cells_collection[0].number_of_contacts == 0);
    reset_forces(cells);
    calculate_force_due_contacts(cells, radius, friction_coefficient);
    compute_new_position(cells, timestep, gravity,1);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * 3.0 = 26.5
    CHECK(cells.cells_collection[0].x == Approx(34.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * 3.0 + (-9.81) * 3.0 * 3.0 * 0.5 = -4.905e-08
    CHECK(cells.cells_collection[0].y == Approx(-19.144999999999996).margin(1e-10));
    // same velocity
    CHECK(cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 3
    CHECK(cells.cells_collection[0].vy == Approx(-29.43).margin(1e-8));
}

TEST_CASE("Two timesteps of 2 living cell and a wall cell")
{

    // Create a 200 cells randomly positioned
    cells_type cells;
    grid_type grid;

    double radius = 4.24;

    double gravity = -9.81;

    cell cell_;
    double width = 100.0;
    double height = 100.0;

    grid.grid_initialize(radius,width,height);

    cell_.x = 25.0;
    cell_.y = 25.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.k = 1.0e6;
    cell_.damp_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.is_wall = -1;
    cells.cells_collection.push_back(cell_);
    cells.number_of_cells++;

    cell_.x = 31.0;
    cell_.y = 31.0;
    cell_.vx = 0.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.k = 2.0e6;
    cell_.damp_ratio = 0.2;
    cell_.mass = 2.0;
    cell_.is_wall = 0;
    cells.cells_collection.push_back(cell_);
    cells.number_of_cells++;

    add_cells_to_grid(cells, grid);

    double timestep = .3;
    double friction_coefficient = 0.0;

    find_contacts_grid(cells,grid, radius, false,width,height);
    REQUIRE(cells.cells_collection[0].number_of_contacts == 0);
    reset_forces(cells);
    calculate_force_due_contacts(cells, radius, friction_coefficient);
    // forces must be zero at this point
    CHECK(cells.cells_collection[0].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.cells_collection[0].fy[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.cells_collection[1].fx[0]== Approx(0.0).margin(1e-10));
    CHECK(cells.cells_collection[1].fy[0] == Approx(0.0).margin(1e-10));
    compute_new_position(cells, timestep, gravity,1);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * .3 = 26.5
    CHECK(cells.cells_collection[0].x == Approx(25.9).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * .3 + (-9.81) * .3 * .3 * 0.5 = -4.905e-08
    CHECK(cells.cells_collection[0].y == Approx(24.55855).margin(1e-10));
    // same velocity
    CHECK(cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 1e-4
    CHECK(cells.cells_collection[0].vy == Approx(-2.943).margin(1e-8));

    // same, since the cell is a wall cell
    CHECK(cells.cells_collection[1].x == Approx(31.0).margin(1e-8));
    CHECK(cells.cells_collection[1].y == Approx(31.0).margin(1e-8));
    CHECK(cells.cells_collection[1].vx == Approx(0.0).margin(1e-8));
    CHECK(cells.cells_collection[1].vy == Approx(0.0).margin(1e-8));

    // NEW TIMESTEP
    find_contacts_grid(cells,grid, radius, false,width,height);
    // WE MUST HAVE ONE CONTACT HERE
    REQUIRE(cells.cells_collection[0].number_of_contacts == 1);
    reset_forces(cells);
    calculate_force_due_contacts(cells, radius, friction_coefficient);

    CHECK(cells.cells_collection[0].fx[0] == Approx(-245579.1171289711).margin(1e-7));
    CHECK(cells.cells_collection[0].fy[0] == Approx(-310173.6478491002).margin(1e-7));
    CHECK(cells.cells_collection[1].fx[0] == Approx(245579.1171289711).margin(1e-7));
    CHECK(cells.cells_collection[1].fy[0] == Approx(310173.6478491002).margin(1e-7));

    compute_new_position(cells, timestep, gravity,1);

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 25.9 + 3.0 * .3 + (-248614.96935456092)/1.0 *.3*.3*0.5= 26.5
    CHECK(cells.cells_collection[0].x == Approx(-11024.2602708037).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 24.55855 + -2.943 * .3 + (-9.81 -307770.0203505543) * .3 * .3 * 0.5 = -4.905e-08
    CHECK(cells.cells_collection[0].y == Approx(-13934.5799532095).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 3.0 + (-248614.96935456092) * .3
    CHECK(cells.cells_collection[0].vx == Approx(-73670.7351386913).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -2.943 + (-9.81 -307770.0203505543) * .3
    CHECK(cells.cells_collection[0].vy == Approx(-93057.9803547301).margin(1e-8));

    // same, since the cell is a wall cell
    CHECK(cells.cells_collection[1].x == Approx(31.0).margin(1e-8));
    CHECK(cells.cells_collection[1].y == Approx(31.0).margin(1e-8));
    CHECK(cells.cells_collection[1].vx == Approx(0.0).margin(1e-8));
    CHECK(cells.cells_collection[1].vy == Approx(0.0).margin(1e-8));
}


TEST_CASE("TWO LIVING CELLS COLIDING WITH LIMITS")
{
    double radius = 3;

    int number_of_cells_alive = 1;
    double output_interval = 1;

    double width = 40.0;
    double height = 40.0;
    double friction_coefficient = 0.3;

    Environment env = Environment(radius, width, height,2.0,0.001, friction_coefficient,true);

    env.initialize_living_cells(number_of_cells_alive, 1e4, 1.1e4, 0.1, 0.1, 10.0, 200.0, width, height,1);

    env.cells.cells_collection[0].x = 9.0;
    env.cells.cells_collection[0].y = 15.0;
    env.cells.cells_collection[0].vx = -9.0;
    env.cells.cells_collection[0].vy = -8.0;
    env.cells.cells_collection[0].k = 10193.304239020967;
    env.cells.cells_collection[0].mass = 1.0;
    env.cells.cells_collection[0].damp_ratio = 0.10000000000000001;

    env.cells.cells_collection[1].x = 25.0;
    env.cells.cells_collection[1].y = 31.0;
    env.cells.cells_collection[1].vx = 9.0;
    env.cells.cells_collection[1].vy = 8.0;
    env.cells.cells_collection[1].k = 10193.304239020967;
    env.cells.cells_collection[1].mass = 1.0;
    env.cells.cells_collection[1].damp_ratio = 0.10000000000000001;

    env.start_simulation(0.1, 1e-3,2.0,1,false);

    CHECK(env.cells.cells_collection[0].x == Approx(4.9878999118).margin(1e-5));

    CHECK(env.cells.cells_collection[0].y == Approx(6.8257076293).margin(1e-5));

    CHECK(env.cells.cells_collection[0].vx == Approx(0.008725002).margin(1e-5));

    CHECK(env.cells.cells_collection[0].vy == Approx(-0.4888244935).margin(1e-5));

    CHECK(env.cells.cells_collection[1].x == Approx(25.0).margin(1e-5));
    CHECK(env.cells.cells_collection[1].y == Approx(31.0).margin(1e-5));
    CHECK(env.cells.cells_collection[1].vx == Approx(9.0).margin(1e-5));
    CHECK(env.cells.cells_collection[1].vy == Approx(8.0).margin(1e-5));
}
