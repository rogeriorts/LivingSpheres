#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../include/LivingSpheres/simulator.h"
using namespace std;

TEST_CASE("Two timesteps of 2 living cells")
{

    // Create a 200 cells randomly positioned
    cells_type cells;
    cell_contacts contacts;

    double radius = 4.24;

    double gravity = -9.81;

    cells.x.push_back(25.0);
    cells.y.push_back(25.0);
    cells.vx.push_back(3.0);
    cells.vy.push_back(0.0);
    cells.fx.push_back(0.0);
    cells.fy.push_back(0.0);
    cells.k.push_back(1.0e6);
    cells.damp_ratio.push_back(0.1);
    cells.mass.push_back(1.0);
    cells.is_wall.push_back(-1);
    cells.number_of_cells++;

    cells.x.push_back(31.0);
    cells.y.push_back(31.0);
    cells.vx.push_back(0.0);
    cells.vy.push_back(-3.0);
    cells.fx.push_back(0.0);
    cells.fy.push_back(0.0);
    cells.k.push_back(2.0e6);
    cells.damp_ratio.push_back(0.2);
    cells.mass.push_back(2.0);
    cells.is_wall.push_back(-1);
    cells.number_of_cells++;

    double timestep = .3;
    double friction_coefficient = 0.3;

    find_contacts(cells, contacts, radius, false);
    REQUIRE(contacts.number_of_contacts == 0);
    reset_forces(cells);
    calculate_force_due_contacts(cells, contacts, radius, friction_coefficient);
    // forces must be zero at this point
    CHECK(cells.fx[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.fy[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.fx[1] == Approx(0.0).margin(1e-10));
    CHECK(cells.fy[1] == Approx(0.0).margin(1e-10));
    compute_new_position(cells, timestep, gravity);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * .3 = 26.5
    CHECK(cells.x[0] == Approx(25.9).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * .3 + (-9.81) * .3 * .3 * 0.5
    CHECK(cells.y[0] == Approx(24.55855).margin(1e-10));
    // same velocity
    CHECK(cells.vx[0] == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 1e-4
    CHECK(cells.vy[0] == Approx(-2.943).margin(1e-8));

    // same, since there is no vx at this point
    CHECK(cells.x[1] == Approx(31.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y =  31.0 + (-3) * .3 + (-9.81) * .3 * .3 * 0.5
    CHECK(cells.y[1] == Approx(29.65855).margin(1e-8));
    CHECK(cells.vx[1] == Approx(0.0).margin(1e-8));
    // same velocity since there is no acceleration
    // vy = vy_0 + g * t
    // vy = -3 + -9.81 * .3
    CHECK(cells.vy[1] == Approx(-5.943).margin(1e-8));

    // NEW TIMESTEP
    find_contacts(cells, contacts, radius, false);
    // WE MUST HAVE ONE CONTACT HERE
    REQUIRE(contacts.number_of_contacts == 1);
    reset_forces(cells);
    calculate_force_due_contacts(cells, contacts, radius, friction_coefficient);

    CHECK(cells.fx[0] == Approx(-1347178.0769591301).margin(1e-7));
    CHECK(cells.fy[0] == Approx(-1347178.0769591301).margin(1e-7));
    CHECK(cells.fx[1] == Approx(1347178.0769591301).margin(1e-7));
    CHECK(cells.fy[1] == Approx(1347178.0769591301).margin(1e-7));

    compute_new_position(cells, timestep, gravity);

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 25.9 + 3.0 * .3 + (-1347178.0769591301)/1.0 *.3*.3*0.5= 26.5
    CHECK(cells.x[0] == Approx(-60596.21346316085).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 24.55855 + -2.943 * .3 + (-9.81 -1347178.0769591301) * .3 * .3 * 0.5
    CHECK(cells.y[0] == Approx(-60599.77926316086).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 3.0 + (-1347178.0769591301) * .3
    CHECK(cells.vx[0] == Approx(-404150.42308773904).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -2.943 + (-9.81 -1347178.0769591301) * .3
    CHECK(cells.vy[0] == Approx(-404159.3090877391).margin(1e-8));

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 31.0 + 0.0 * .3 + (1347178.0769591301)/2.0 *.3*.3*0.5
    CHECK(cells.x[1] == Approx(30342.506731580426).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 29.65855 + -5.943 * .3 + (-9.81 +1347178.0769591301/2.0) * .3 * .3 * 0.5
    CHECK(cells.y[1] == Approx(30338.940931580422).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 0.0 + (1347178.0769591301/2.0) * .3
    CHECK(cells.vx[1] == Approx(202076.71154386952).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -5.943 + (-9.81 +1347178.0769591301)/2 * .3
    CHECK(cells.vy[1] == Approx(202069.2970438695).margin(1e-8));
}

TEST_CASE("One cell falling due gravity")
{

    // Create a 200 cells randomly positioned
    cells_type cells;
    cell_contacts contacts;

    double radius = 4.24;

    double gravity = -9.81;

    cells.x.push_back(25.0);
    cells.y.push_back(25.0);
    cells.vx.push_back(3.0);
    cells.vy.push_back(0.0);
    cells.fx.push_back(0.0);
    cells.fy.push_back(0.0);
    cells.k.push_back(1.0e6);
    cells.damp_ratio.push_back(0.1);
    cells.mass.push_back(1.0);
    cells.is_wall.push_back(-1);
    cells.number_of_cells++;

    double timestep = 3.0;
    double friction_coefficient = 0.3;

    find_contacts(cells, contacts, radius, false);
    REQUIRE(contacts.number_of_contacts == 0);
    reset_forces(cells);
    calculate_force_due_contacts(cells, contacts, radius, friction_coefficient);
    compute_new_position(cells, timestep, gravity);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * 3.0 = 26.5
    CHECK(cells.x[0] == Approx(34.0).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * 3.0 + (-9.81) * 3.0 * 3.0 * 0.5 = -4.905e-08
    CHECK(cells.y[0] == Approx(-19.144999999999996).margin(1e-10));
    // same velocity
    CHECK(cells.vx[0] == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 3
    CHECK(cells.vy[0] == Approx(-29.43).margin(1e-8));
}

TEST_CASE("Two timesteps of 2 living cell and a wall cell")
{

    // Create a 200 cells randomly positioned
    cells_type cells;
    cell_contacts contacts;

    double radius = 4.24;

    double gravity = -9.81;

    cells.x.push_back(25.0);
    cells.y.push_back(25.0);
    cells.vx.push_back(3.0);
    cells.vy.push_back(0.0);
    cells.fx.push_back(0.0);
    cells.fy.push_back(0.0);
    cells.k.push_back(1.0e6);
    cells.damp_ratio.push_back(0.1);
    cells.mass.push_back(1.0);
    cells.is_wall.push_back(-1);
    cells.number_of_cells++;

    cells.x.push_back(31.0);
    cells.y.push_back(31.0);
    cells.vx.push_back(0.0);
    cells.vy.push_back(0.0);
    cells.fx.push_back(0.0);
    cells.fy.push_back(0.0);
    cells.k.push_back(2.0e6);
    cells.damp_ratio.push_back(0.2);
    cells.mass.push_back(2.0);
    cells.is_wall.push_back(0);
    cells.number_of_cells++;

    double timestep = .3;
    double friction_coefficient = 0.0;

    find_contacts(cells, contacts, radius, false);
    REQUIRE(contacts.number_of_contacts == 0);
    reset_forces(cells);
    calculate_force_due_contacts(cells, contacts, radius, friction_coefficient);
    // forces must be zero at this point
    CHECK(cells.fx[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.fy[0] == Approx(0.0).margin(1e-10));
    CHECK(cells.fx[1] == Approx(0.0).margin(1e-10));
    CHECK(cells.fy[1] == Approx(0.0).margin(1e-10));
    compute_new_position(cells, timestep, gravity);

    // x = x0 + vx * t
    // x = 25.0 + 3.0 * .3 = 26.5
    CHECK(cells.x[0] == Approx(25.9).margin(1e-8));
    // y = y0 + vy_0 * t + g * t * t * 0.5;
    // y = 25.0 + 0.0 * .3 + (-9.81) * .3 * .3 * 0.5 = -4.905e-08
    CHECK(cells.y[0] == Approx(24.55855).margin(1e-10));
    // same velocity
    CHECK(cells.vx[0] == Approx(3.0).margin(1e-8));
    // vy = vy_0 + g * t
    // vy = 0.0 + -9.81 * 1e-4
    CHECK(cells.vy[0] == Approx(-2.943).margin(1e-8));

    // same, since the cell is a wall cell
    CHECK(cells.x[1] == Approx(31.0).margin(1e-8));
    CHECK(cells.y[1] == Approx(31.0).margin(1e-8));
    CHECK(cells.vx[1] == Approx(0.0).margin(1e-8));
    CHECK(cells.vy[1] == Approx(0.0).margin(1e-8));

    // NEW TIMESTEP
    find_contacts(cells, contacts, radius, false);
    // WE MUST HAVE ONE CONTACT HERE
    REQUIRE(contacts.number_of_contacts == 1);
    reset_forces(cells);
    calculate_force_due_contacts(cells, contacts, radius, friction_coefficient);

    CHECK(cells.fx[0] == Approx(-245579.1171289711).margin(1e-7));
    CHECK(cells.fy[0] == Approx(-310173.6478491002).margin(1e-7));
    CHECK(cells.fx[1] == Approx(245579.1171289711).margin(1e-7));
    CHECK(cells.fy[1] == Approx(310173.6478491002).margin(1e-7));

    compute_new_position(cells, timestep, gravity);

    // x = x0 + vx_0 * t + fx/m * t^2 * 0.5
    // x = 25.9 + 3.0 * .3 + (-248614.96935456092)/1.0 *.3*.3*0.5= 26.5
    CHECK(cells.x[0] == Approx(-11024.2602708037).margin(1e-5));
    // y = y0 + vy_0 * t + (g + fy/m) * t * t * 0.5;
    // y = 24.55855 + -2.943 * .3 + (-9.81 -307770.0203505543) * .3 * .3 * 0.5 = -4.905e-08
    CHECK(cells.y[0] == Approx(-13934.5799532095).margin(1e-10));
    // vx = vx_0 + (fx/m) * t
    // vx = 3.0 + (-248614.96935456092) * .3
    CHECK(cells.vx[0] == Approx(-73670.7351386913).margin(1e-8));
    // vy = vy_0 + (g + fy/m) * t
    // vy = -2.943 + (-9.81 -307770.0203505543) * .3
    CHECK(cells.vy[0] == Approx(-93057.9803547301).margin(1e-8));

    // same, since the cell is a wall cell
    CHECK(cells.x[1] == Approx(31.0).margin(1e-8));
    CHECK(cells.y[1] == Approx(31.0).margin(1e-8));
    CHECK(cells.vx[1] == Approx(0.0).margin(1e-8));
    CHECK(cells.vy[1] == Approx(0.0).margin(1e-8));
}

TEST_CASE("TWO LIVING CELLS COLIDING WITH LIMITS")
{
    double radius = 3;

    int number_of_cells_alive = 2;
    double output_interval = 1;

    double width = 40.0;
    double height = 40.0;
    double friction_coefficient = 0.3;

    Environment env = Environment(radius, width, height, friction_coefficient,true);

    env.initialize_living_cells(number_of_cells_alive, 1e4, 1.1e4, 0.1, 0.1, 10.0, 200.0, width, height);

    env.cells.x[0] = 9.0;
    env.cells.y[0] = 15.0;
    env.cells.vx[0] = -9.0;
    env.cells.vy[0] = -8.0;
    env.cells.k[0] = 10193.304239020967;
    env.cells.mass[0] = 1.0;
    env.cells.damp_ratio[0] = 0.10000000000000001;

    env.cells.x[1] = 25.0;
    env.cells.y[1] = 31.0;
    env.cells.vx[1] = 9.0;
    env.cells.vy[1] = 8.0;
    env.cells.k[1] = 10193.304239020967;
    env.cells.mass[1] = 1.0;
    env.cells.damp_ratio[1] = 0.10000000000000001;

    env.start_simulation(0.1, 2,false);

    CHECK(env.cells.x[0] == Approx(4.8707516304).margin(1e-5));

    CHECK(env.cells.y[0] == Approx(6.3189700679).margin(1e-5));

    CHECK(env.cells.vx[0] == Approx(-0.0037469552).margin(1e-5));

    CHECK(env.cells.vy[0] == Approx(-1.0324224535).margin(1e-5));

    CHECK(env.cells.x[1] == Approx(33.6367230952).margin(1e-5));
    CHECK(env.cells.y[1] == Approx(30.1710589546).margin(1e-5));
    CHECK(env.cells.vx[1] == Approx(-5.2915778695).margin(1e-5));
    CHECK(env.cells.vy[1] == Approx(-7.3355395629).margin(1e-5));
}