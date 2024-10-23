/******************************************************************************

     MIT License

    Copyright (c) 2024 Francesco Zattoni

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

*********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "sequential/functions.h"
#include "sequential/physics.h"
#include "sequential/simulation_utils.h"

#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

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
    
    struct timeval t_start_1, t_end_1, t_start_2, t_end_2;

    /**
     * Placing the cells and antigens in the grid and placing only free cells in the next iteration grid.
     */
    generation(&grid, options);
    free_grid(&next_grid);

    //Save cells at the start in a file
    save_grid(&grid, "start.csv");

    gettimeofday(&t_start_1, NULL);
    simulation(&grid, &next_grid, options);
    gettimeofday(&t_end_1, NULL);

    //Save cells in the middle of the simulation, before inserting new antigens.
    save_grid(&grid, "mid.csv");
    insert_antigens(&grid);

    gettimeofday(&t_start_2, NULL);
    simulation(&grid, &next_grid, options);
    gettimeofday(&t_end_2, NULL);

    printf("T: %lf\n", (WALLTIME(t_end_1) - WALLTIME(t_start_1)) + (WALLTIME(t_end_2) - WALLTIME(t_start_2)));

    //Save cells at the end in a file
    save_grid(&grid, "end.csv");

    free(grid.matrix);
    free(next_grid.matrix);
    return 0;
}