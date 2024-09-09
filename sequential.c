#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "functions.h"
#include "physics.h"
#include "simulation_utils.h"

int main(int argc, char const *argv[])
{
    srand(time(NULL));

    Options options;
    read_parameters(&options, argv, argc - 1);

    Cell* grid = (Cell*)malloc(sizeof(Cell) * options.grid_size * options.grid_size),
        *next_grid = (Cell*)malloc(sizeof(Cell) * options.grid_size * options.grid_size),
        *temp;
    generation(grid, options);

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
    }

    free(grid);
    free(next_grid);
    return 0;
}