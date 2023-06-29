#include <vector>

struct cell
{
    double x;
    double y;
    double vx;
    double vy;
    std::vector<double> fx;
    std::vector<double> fy;
    std::vector<double> transferred_heat;


    double mass;
    double volume;
    double damping_ratio;
    double temperature;
    
    //material properties
    double spring_coefficient;
    double density;
    double heat_conduction_coefficient;
    double specific_heat;

    bool fixed = false;
    
    int is_wall;
    std::vector<int> contacts;
    std::vector<double> contact_vector_x;
    std::vector<double> contact_vector_y;
    int number_of_contacts = 0;
    int grid_block = -1;

    void add_contact(int cell_b_, double contact_vector_x_,
                     double contact_vector_y_)
    {
        number_of_contacts++;

        if (number_of_contacts - 1 >= contacts.size())
        {
            contacts.push_back(cell_b_);
            contact_vector_x.push_back(contact_vector_x_);
            contact_vector_y.push_back(contact_vector_y_);
        }
        else
        {
            contacts[number_of_contacts - 1] = cell_b_;
            contact_vector_x[number_of_contacts - 1] = contact_vector_x_;
            contact_vector_y[number_of_contacts - 1] = contact_vector_y_;
        }
    }
};

struct cells_type
{

    std::vector<cell> cells_collection;
    int number_of_cells = 0;
};