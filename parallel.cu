#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>

#include "parallel/functions.h"
#include "parallel/physics.h"
#include "parallel/simulation_utils.h"

/**
 * Individual timestep ran by each thread, it will be the main kernel of the simulation.
 * @param grid grid of the current iteration
 * @param next_grid grid of the next iteration
 */
__global__ void timestep(Grid* grid, Grid* next_grid, curandState* rand_state)
{
    /**
     * Selecting a cell in position (x, y), executing its action and making it move.
     */
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < grid->size  && y < grid->size)
    {
        Vector position = {(float) x,(float) y};
        Cell* cell = device_access_grid(grid, position);
        cell->action(cell, grid, next_grid);  
        //movement(cell, next_grid, rand_state);
    }
}

/**
 * Actual simulation function that for n timesteps updates cells in the grid.
 * @param grid grid of the current iteration
 * @param next_grid grid of the next iteration
 * @param options simulation options
 */
__host__ void simulation(Grid* grid, Grid* next_grid, Options options)
{
    dim3 block_size = {32, 32, 1};
    unsigned int grid_length_x = (grid->size) / block_size.x + 1;
    unsigned int grid_length_y = (grid->size) / block_size.y + 1;
    dim3 number_blocks = {grid_length_x, grid_length_y, 1};
    dim3 thread_per_block = {block_size.x, block_size.y, 1};

    curandState *rand_state;
    cudaMalloc(&rand_state, sizeof(curandState));

    for (int t = 0; t < TIMESTEPS; t++)
    {
        timestep<<<number_blocks, thread_per_block>>>(grid, next_grid, rand_state);
        /**
         * Swapping grids, so that for every iteration it is used the same grid variable.
         */
        swap_grids<<<number_blocks, thread_per_block>>>(grid, next_grid, options.grid_size);
    }

    cudaFree(rand_state);
}

int main(int argc, char const *argv[])
{
    srand(time(NULL));

    /**
     * Initialization of simulation options and allocation of grids' memory.
     */
    Options options;
    read_parameters(&options, argv, argc - 1);

    Grid h_grid, h_next_grid,
        d_grid, d_next_grid;

    int malloc_size = sizeof(Cell) * options.grid_size * options.grid_size;

    d_grid.size = d_next_grid.size = h_grid.size = h_next_grid.size = options.grid_size;
    h_grid.matrix = (Cell*)malloc(malloc_size),
    h_next_grid.matrix = (Cell*)malloc(malloc_size);
    cudaMalloc((void**) &(d_grid.matrix), malloc_size);
    cudaMalloc((void**) &(d_next_grid.matrix), malloc_size);

    /**
     * Placing the cells and antigens in the grid and placing only free cells in the next iteration grid.
     */
    generation(&h_grid, options);
    free_grid(&h_next_grid);
    cudaMemcpy((void*) d_grid.matrix, (void*) h_grid.matrix, malloc_size, cudaMemcpyHostToDevice);
    cudaMemcpy((void*) d_next_grid.matrix, (void*) h_next_grid.matrix, malloc_size, cudaMemcpyHostToDevice);

    //Save cells at the start in a file
    const char* start_file = "grids/start.csv"; 
    save_grid(&h_grid, start_file);

    simulation(&d_grid, &d_next_grid, options);

    //Save cells in the middle of the simulation, before inserting new antigens.
    cudaMemcpy((void*) h_grid.matrix, (void*) d_grid.matrix, malloc_size, cudaMemcpyDeviceToHost);
    const char* mid_file = "grids/mid.csv";
    save_grid(&h_grid, mid_file);
    insert_antigens(&h_grid);
    cudaMemcpy((void*) d_grid.matrix, (void*) h_grid.matrix, malloc_size, cudaMemcpyHostToDevice);

    simulation(&d_grid, &d_next_grid, options);

    //Save cells at the end in a file
    cudaMemcpy((void*) h_grid.matrix, (void*) d_grid.matrix, malloc_size, cudaMemcpyDeviceToHost);
    const char* end_file = "grids/end.csv";
    save_grid(&h_grid, end_file);

    free(h_grid.matrix);
    free(h_next_grid.matrix);
    cudaFree(d_grid.matrix);
    cudaFree(d_next_grid.matrix);
    return 0;
}