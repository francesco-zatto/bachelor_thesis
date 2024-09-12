#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "sequential/functions.h"
#include "sequential/physics.h"
#include "sequential/simulation_utils.h"

/**
 * Actual simulation function that for n timesteps updates cells in the grid.
 * @param grid grid of the current iteration
 * @param next_grid grid of the next iteration
 * @param options simulation options
 */
void simulation(Grid* grid, Grid* next_grid, Options options)
{
    for (int t = 0; t < TIMESTEPS; t++)
    {
        for (int i = 0; i < options.grid_size; i++)
        {
            for (int j = 0; j < options.grid_size; j++)
            {
                /**
                 * Selecting a cell in position (i, j), executing its action and making it move.
                 */
                Vector position = {i, j};
                Cell* cell = access_grid(grid, position);
                cell->action(cell, grid, next_grid);  
                movement(cell, next_grid);
            }
        }
        /**
         * Swapping grids, so that for every iteration it is used the same grid variable.
         */
        swap_grids(grid, next_grid, options.grid_size);
    }
}

int main(int argc, char const *argv[])
{
    srand(time(NULL));

    /**
     * Initialization of simulation options and allocation of grids' memory.
     */
    Options options;
    read_parameters(&options, argv, argc - 1);

    Grid grid, next_grid;
    grid.size = next_grid.size = options.grid_size;
    grid.matrix = (Cell*)malloc(sizeof(Cell) * options.grid_size * options.grid_size),
    next_grid.matrix = (Cell*)malloc(sizeof(Cell) * options.grid_size * options.grid_size);
    
    /**
     * Placing the cells and antigens in the grid and placing only free cells in the next iteration grid.
     */
    generation(&grid, options);
    free_grid(&next_grid);

    //Save cells at the start in a file
    save_grid(&grid, "grids/start.csv");

    simulation(&grid, &next_grid, options);

    //Save cells in the middle of the simulation, before inserting new antigens.
    save_grid(&grid, "grids/mid.csv");
    insert_antigens(&grid);

    simulation(&grid, &next_grid, options);

    //Save cells at the end in a file
    save_grid(&grid, "grids/end.csv");

    free(grid.matrix);
    free(next_grid.matrix);
    return 0;
}