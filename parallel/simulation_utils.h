#ifndef _SIM_UTILS_HEADER
#define _SIM_UTILS_HEADER

#include "../cell.h"

/**
 * 2D Vector with coordinates set to 0.
 */
extern const Vector NULL_VECTOR;

/**
 * Cell with type set to free. Used to free the next grid.
 */
extern const Cell FREE_CELL;

/**
 * Timesteps of the HIS simulation.
 */
extern const int TIMESTEPS;

/**
 * Size of grid size.
 */
extern const int GRID_SIZE;

/**
 * Default number of B lymphocytes of the simulation.
 */
extern const int CELLS_B_NUMBER;

/**
 * Default number of T helper lymphocytes of the simulation.
 */
extern const int CELLS_T_NUMBER;

/**
 * Default number of antigens of the simulation.
 */
extern const int AG_NUMBER;

/**
 * Default number of new inserted antigens in the middle of the simulation.
 */
extern const int NEW_AG_NUMBER;

/**
 * Function to read command line parameters to change default values,
 * like numbers of cells or grid size.
 * @param options output with the setted options
 * @param parameters command line parameters
 * @param n number of parameters to save
 */
__host__ void read_parameters(Options* options, const char* parameters[], int n);

/**
 * Function to generate the grid, taking generation's options.
 * @param grid grid to generate
 * @param options generation options, like B cells quantity
 */
__host__ void generation(Grid* grid, Options options);

/**
 * Function to create a single cell given its position and its type.
 * @param cell output cell created by the function
 * @param position vector with grid position of the new cell
 * @param type cell type
 */
__host__ void create_cell(Cell* cell, Vector position, Type type);

/**
 * Function to extract the type of the current cell, depending on the options.
 * @param options generation options
 * @return type of the current cell
 */
__host__ Type extract_type(Options options);

/**
 * Function to swap old grid and the new one to continue the simulation with the same grid.
 * @param old_grid old grid of the current timestep
 * @param new_grid new grid of the next timestep
 * @param size grid's side length
 */
__global__ void swap_grids(Grid* old_grid, Grid* new_grid, int size);

/**
 * Function to fill a grid with only free cells. It does not deallocate grid memory.
 * @param grid output grid that will contain only free cells
 */
__host__ void free_grid(Grid* grid);

/**
 * Function to save every information about the cells inside the given grid in the current state
 * @param grid grid to save in external file
 * @param filename name of the file
 */
__host__ void save_grid(Grid* grid, const char* filename);

/**
 * Function to insert new antigens in the simulation's grid to test HIS memory.
 * @param grid grid which antigens will be inserted in
 */
__host__ void insert_antigens(Grid* grid);

#endif