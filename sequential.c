#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "functions.h"
#include "physics.h"
#include "simulation_utils.h"

void simulation(Grid* grid, Grid* next_grid, Options options)
{
    for (int t = 0; t < TIMESTEPS; t++)
    {
        for (int i = 0; i < options.grid_size; i++)
        {
            for (int j = 0; j < options.grid_size; j++)
            {
                Vector position = {i, j};
                Cell* cell = access_grid(grid, position);
                cell->action(cell, grid, next_grid);  
                movement(cell, next_grid);
            }
        }
        swap_grids(grid, next_grid, options.grid_size);
        if (t == TIMESTEPS / 2)
        {
            save_grid(grid, "grids/mid.csv");
        }
    }
}

int main(int argc, char const *argv[])
{
    srand(time(NULL));

    Options options;
    read_parameters(&options, argv, argc - 1);

    Grid grid, next_grid;
    grid.size = next_grid.size = options.grid_size;
    grid.matrix = (Cell*)malloc(sizeof(Cell) * options.grid_size * options.grid_size),
    next_grid.matrix = (Cell*)malloc(sizeof(Cell) * options.grid_size * options.grid_size);
    
    generation(&grid, options);
    free_grid(&next_grid);

    save_grid(&grid, "grids/start.csv");

    simulation(&grid, &next_grid, options);

    save_grid(&grid, "grids/end.csv");

    free(grid.matrix);
    free(next_grid.matrix);
    return 0;
}