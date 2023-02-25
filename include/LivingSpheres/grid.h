#include <vector>
#include <cmath>

struct grid_block_type
{
    std::vector<int> cell_ids; 
    int next_position_on_block = 0;

};

struct grid_type{
    std::vector<grid_block_type> blocks;
    double block_size_x,block_size_y;
    int nx,ny;
    
    void grid_initialize(double radius, double width, double height,double grid_multiplier=3.0)
    {
        block_size_x = grid_multiplier*radius;
        block_size_y = block_size_x;

        nx = (int)std::floor(width / block_size_x);
        block_size_x = width / (double)nx;
        
        ny = (int)std::floor(height / block_size_y);
        block_size_y = height / (double)ny;

        blocks.resize(nx*ny);

    }

    bool add_cell_to_block(double cell_x, double cell_y, int cell_index)
    {
        if (cell_x > block_size_x * (double)nx || cell_y > block_size_y * (double)ny
        ||cell_x < 0.0 || cell_y < 0.0) 
            return false;
        
        int i = (int) std::floor(cell_x / block_size_x);
        int j = (int) std::floor(cell_y / block_size_y);
        
        int block_index = i + nx * j;

        if(blocks[block_index].next_position_on_block >= blocks[block_index].cell_ids.size())
        {
            blocks[block_index].cell_ids.push_back(cell_index);
        }
        else
        {
            blocks[block_index].cell_ids[blocks[block_index].next_position_on_block] = cell_index;
        }
        blocks[block_index].next_position_on_block++;

        return true;
    }
    
    void reset_grid()
    {
        for(int i = 0; i < blocks.size();i++)
            blocks[i].next_position_on_block = 0;
    }

};