#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../include/LivingSpheres/simulator.h"
using namespace std;

TEST_CASE("Two timesteps of 2 living cells")
{

    // Create a 200 cells randomly positioned

    double radius = 4.24;

    double gravity = -9.81;

    cell cell_;
    double width = 100.0;
    double height = 100.0;

    double timestep = .3;
    double friction_coefficient = 0.3;

    Simulator sim = Simulator(
        radius, width, height, 3.0,
        0.01, friction_coefficient, 1, false, gravity);

    sim.timestep = timestep;
    cell_.x = 25.0;
    cell_.y = 25.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.spring_coefficient = 1.0e6;
    cell_.damping_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.volume = (4.0 / 3.0 * M_PI * radius * radius * radius);
    cell_.density = cell_.mass / cell_.volume;
    cell_.is_wall = -1;
    sim.cells.cells_collection.push_back(cell_);
    sim.cells.number_of_cells++;

    cell_.x = 31.0;
    cell_.y = 31.0;
    cell_.vx = 0.0;
    cell_.vy = -3.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.spring_coefficient = 2.0e6;
    cell_.damping_ratio = 0.2;
    cell_.mass = 2.0;
    cell_.is_wall = -1;
    sim.cells.cells_collection.push_back(cell_);
    sim.cells.number_of_cells++;

    sim.add_cells_to_grid();

    sim.find_contacts_grid();

    sim.reset_forces();
    sim.calculate_force_due_contacts();
    // forces must be zero at this point
    CHECK(sim.cells.cells_collection[0].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(sim.cells.cells_collection[0].fy[0] == Approx(0.0).margin(1e-10));
    CHECK(sim.cells.cells_collection[1].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(sim.cells.cells_collection[1].fy[0] == Approx(0.0).margin(1e-10));
    sim.compute_new_position();

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * .3 = 26.5
    CHECK(sim.cells.cells_collection[0].x == Approx(25.9).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * .3 + (-9.81) * .3 * .3 * 0.5
    CHECK(sim.cells.cells_collection[0].y == Approx(24.55855).margin(1e-10));
    // same velocity
    CHECK(sim.cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 1e-4
    CHECK(sim.cells.cells_collection[0].vy == Approx(-2.943).margin(1e-8));

    // same, since there is no vx at this point
    CHECK(sim.cells.cells_collection[1].x == Approx(31.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y =  31.0 + (-3) * .3 + (-9.81) * .3 * .3 * 0.5
    CHECK(sim.cells.cells_collection[1].y == Approx(29.65855).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].vx == Approx(0.0).margin(1e-8));
    // same velocity since there is no acceleration
    // vy = vy_0 + g * t
    // vy = -3 + -9.81 * .3
    CHECK(sim.cells.cells_collection[1].vy == Approx(-5.943).margin(1e-8));

    // NEW TIMESTEP
    sim.find_contacts_grid();
    // WE MUST HAVE ONE CONTACT HERE
    REQUIRE(sim.cells.cells_collection[0].number_of_contacts == 1);
    sim.reset_forces();
    sim.calculate_force_due_contacts();

    CHECK(sim.cells.cells_collection[0].fx[0] == Approx(-1347178.0769591301).margin(1e-7));
    CHECK(sim.cells.cells_collection[0].fy[0] == Approx(-1347178.0769591301).margin(1e-7));
    CHECK(sim.cells.cells_collection[1].fx[0] == Approx(1347178.0769591301).margin(1e-7));
    CHECK(sim.cells.cells_collection[1].fy[0] == Approx(1347178.0769591301).margin(1e-7));

    sim.compute_new_position();

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 25.9 + 3.0 * .3 + (-1347178.0769591301)/1.0 *.3*.3*0.5= 26.5
    CHECK(sim.cells.cells_collection[0].x == Approx(-60596.21346316085).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 24.55855 + -2.943 * .3 + (-9.81 -1347178.0769591301) * .3 * .3 * 0.5
    CHECK(sim.cells.cells_collection[0].y == Approx(-60599.77926316086).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 3.0 + (-1347178.0769591301) * .3
    CHECK(sim.cells.cells_collection[0].vx == Approx(-404150.42308773904).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -2.943 + (-9.81 -1347178.0769591301) * .3
    CHECK(sim.cells.cells_collection[0].vy == Approx(-404159.3090877391).margin(1e-8));

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 31.0 + 0.0 * .3 + (1347178.0769591301)/2.0 *.3*.3*0.5
    CHECK(sim.cells.cells_collection[1].x == Approx(30342.506731580426).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 29.65855 + -5.943 * .3 + (-9.81 +1347178.0769591301/2.0) * .3 * .3 * 0.5
    CHECK(sim.cells.cells_collection[1].y == Approx(30338.940931580422).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 0.0 + (1347178.0769591301/2.0) * .3
    CHECK(sim.cells.cells_collection[1].vx == Approx(202076.71154386952).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -5.943 + (-9.81 +1347178.0769591301)/2 * .3
    CHECK(sim.cells.cells_collection[1].vy == Approx(202069.2970438695).margin(1e-8));
}

TEST_CASE("One cell falling due gravity")
{

    // Create a 200 cells randomly positioned

    double radius = 4.24;

    double gravity = -9.81;

    cell cell_;
    double width = 100.0;
    double height = 100.0;

    double timestep = 3.0;
    double friction_coefficient = 0.3;

    Simulator sim = Simulator(
        radius, width, height, 3.0,
        0.001, friction_coefficient, 1, false, gravity);
    sim.timestep = timestep;
    cell_.x = 25.0;
    cell_.y = 25.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.spring_coefficient = 1.0e6;
    cell_.damping_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.volume = (4.0 / 3.0 * M_PI * radius * radius * radius);
    cell_.density = cell_.mass / cell_.volume;
    cell_.is_wall = -1;
    sim.cells.cells_collection.push_back(cell_);
    sim.cells.number_of_cells++;

    sim.add_cells_to_grid();

    sim.find_contacts_grid();
    REQUIRE(sim.cells.cells_collection[0].number_of_contacts == 0);
    sim.reset_forces();
    sim.calculate_force_due_contacts();
    sim.compute_new_position();

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * 3.0 = 26.5
    CHECK(sim.cells.cells_collection[0].x == Approx(34.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * 3.0 + (-9.81) * 3.0 * 3.0 * 0.5 = -4.905e-08
    CHECK(sim.cells.cells_collection[0].y == Approx(-19.144999999999996).margin(1e-10));
    // same velocity
    CHECK(sim.cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 3
    CHECK(sim.cells.cells_collection[0].vy == Approx(-29.43).margin(1e-8));
}

TEST_CASE("Two timesteps of 2 living cell and a wall cell")
{

    // Create a 200 cells randomly positioned
    double radius = 4.24;

    double gravity = -9.81;

    cell cell_;
    double width = 100.0;
    double height = 100.0;

    double timestep = .3;
    double friction_coefficient = 0.0;

    Simulator sim = Simulator(
        radius, width, height, 3.0,
        0.001, friction_coefficient, 1, false, gravity);
    
    sim.timestep = timestep;
    
    cell_.x = 25.0;
    cell_.y = 25.0;
    cell_.vx = 3.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.spring_coefficient = 1.0e6;
    cell_.damping_ratio = 0.1;
    cell_.mass = 1.0;
    cell_.volume = (4.0 / 3.0 * M_PI * radius * radius * radius);
    cell_.density = cell_.mass / cell_.volume;
    cell_.is_wall = -1;
    sim.cells.cells_collection.push_back(cell_);
    sim.cells.number_of_cells++;

    cell_.x = 31.0;
    cell_.y = 31.0;
    cell_.vx = 0.0;
    cell_.vy = 0.0;
    cell_.fx.push_back(0.0);
    cell_.fy.push_back(0.0);
    cell_.spring_coefficient = 2.0e6;
    cell_.damping_ratio = 0.2;
    cell_.mass = 2.0;
    cell_.is_wall = 0;
    sim.cells.cells_collection.push_back(cell_);
    sim.cells.number_of_cells++;

    sim.add_cells_to_grid();

    sim.find_contacts_grid();
    REQUIRE(sim.cells.cells_collection[0].number_of_contacts == 0);
    sim.reset_forces();
    sim.calculate_force_due_contacts();
    // forces must be zero at this point
    CHECK(sim.cells.cells_collection[0].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(sim.cells.cells_collection[0].fy[0] == Approx(0.0).margin(1e-10));
    CHECK(sim.cells.cells_collection[1].fx[0] == Approx(0.0).margin(1e-10));
    CHECK(sim.cells.cells_collection[1].fy[0] == Approx(0.0).margin(1e-10));
    sim.compute_new_position();

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * .3 = 26.5
    CHECK(sim.cells.cells_collection[0].x == Approx(25.9).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * .3 + (-9.81) * .3 * .3 * 0.5 = -4.905e-08
    CHECK(sim.cells.cells_collection[0].y == Approx(24.55855).margin(1e-10));
    // same velocity
    CHECK(sim.cells.cells_collection[0].vx == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 1e-4
    CHECK(sim.cells.cells_collection[0].vy == Approx(-2.943).margin(1e-8));

    // same, since the cell is a wall cell
    CHECK(sim.cells.cells_collection[1].x == Approx(31.0).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].y == Approx(31.0).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].vx == Approx(0.0).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].vy == Approx(0.0).margin(1e-8));

    // NEW TIMESTEP
    sim.find_contacts_grid();
    // WE MUST HAVE ONE CONTACT HERE
    REQUIRE(sim.cells.cells_collection[0].number_of_contacts == 1);
    sim.reset_forces();
    sim.calculate_force_due_contacts();

    CHECK(sim.cells.cells_collection[0].fx[0] == Approx(-245579.1171289711).margin(1e-7));
    CHECK(sim.cells.cells_collection[0].fy[0] == Approx(-310173.6478491002).margin(1e-7));
    CHECK(sim.cells.cells_collection[1].fx[0] == Approx(245579.1171289711).margin(1e-7));
    CHECK(sim.cells.cells_collection[1].fy[0] == Approx(310173.6478491002).margin(1e-7));

    sim.compute_new_position();

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 25.9 + 3.0 * .3 + (-248614.96935456092)/1.0 *.3*.3*0.5= 26.5
    CHECK(sim.cells.cells_collection[0].x == Approx(-11024.2602708037).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 24.55855 + -2.943 * .3 + (-9.81 -307770.0203505543) * .3 * .3 * 0.5 = -4.905e-08
    CHECK(sim.cells.cells_collection[0].y == Approx(-13934.5799532095).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 3.0 + (-248614.96935456092) * .3
    CHECK(sim.cells.cells_collection[0].vx == Approx(-73670.7351386913).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -2.943 + (-9.81 -307770.0203505543) * .3
    CHECK(sim.cells.cells_collection[0].vy == Approx(-93057.9803547301).margin(1e-8));

    // same, since the cell is a wall cell
    CHECK(sim.cells.cells_collection[1].x == Approx(31.0).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].y == Approx(31.0).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].vx == Approx(0.0).margin(1e-8));
    CHECK(sim.cells.cells_collection[1].vy == Approx(0.0).margin(1e-8));
}

TEST_CASE("TWO LIVING CELLS COLIDING WITH LIMITS")
{
    double radius = 3;

    int number_of_cells_alive = 2;
    double output_interval = 1;

    double width = 40.0;
    double height = 40.0;
    double friction_coefficient = 0.3;
    double density = 1.0 / (4.0 / 3.0 * M_PI * radius * radius * radius);

    Simulator sim = Simulator(radius, width, height, 2.0, 0.001, friction_coefficient, 1, true);

    sim.initialize_living_cells(number_of_cells_alive, 1e4, 1.1e4, 0.1, 0.1, 10.0, 200.0, width, height,density);

    sim.cells.cells_collection[0].x = 9.0;
    sim.cells.cells_collection[0].y = 15.0;
    sim.cells.cells_collection[0].vx = -9.0;
    sim.cells.cells_collection[0].vy = -8.0;
    sim.cells.cells_collection[0].spring_coefficient = 10193.304239020967;
    sim.cells.cells_collection[0].mass = 1.0;
    sim.cells.cells_collection[0].damping_ratio = 0.10000000000000001;

    sim.cells.cells_collection[1].x = 25.0;
    sim.cells.cells_collection[1].y = 31.0;
    sim.cells.cells_collection[1].vx = 9.0;
    sim.cells.cells_collection[1].vy = 8.0;
    sim.cells.cells_collection[1].spring_coefficient = 10193.304239020967;
    sim.cells.cells_collection[1].mass = 1.0;
    sim.cells.cells_collection[1].damping_ratio = 0.10000000000000001;

    sim.start_simulation(0.1, 1e-3, 2.0, true);

    CHECK(sim.cells.cells_collection[0].x == Approx(4.9847647304 ).margin(1e-5));

    CHECK(sim.cells.cells_collection[0].y == Approx(6.8218351073).margin(1e-5));

    CHECK(sim.cells.cells_collection[0].vx == Approx( 0.0054471666).margin(1e-5));

    CHECK(sim.cells.cells_collection[0].vy == Approx(-0.4931117029 ).margin(1e-5));

    CHECK(sim.cells.cells_collection[1].x == Approx(33.404219618).margin(1e-5));
    CHECK(sim.cells.cells_collection[1].y == Approx(30.2431133273).margin(1e-5));
    CHECK(sim.cells.cells_collection[1].vx == Approx(-5.6497890808).margin(1e-5));
    CHECK(sim.cells.cells_collection[1].vy == Approx(-7.2223913374).margin(1e-5));
}
