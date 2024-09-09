#ifndef _SIM_UTILS_HEADER
#define _SIM_UTILS_HEADER

#include "cell.h"

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
 * Function to read command line parameters to change default values,
 * like numbers of cells or grid size.
 * @param options output with the setted options
 * @param parameters command line parameters
 * @param n number of parameters to save
 */
void read_parameters(Options* options, const char* parameters[], int n);

/**
 * Function to generate the grid, taking generation's options.
 * @param grid grid to generate
 * @param options generation options, like B cells quantity
 */
void generation(Cell* grid, Options options);

/**
 * Function to create a single cell given its position and its type.
 * @param cell output cell created by the function
 * @param position vector with grid position of the new cell
 * @param type cell type
 */
void create_cell(Cell* cell, Vector position, Type type);

/**
 * Function to extract the type of the current cell, depending on the options.
 * @param options generation options
 * @return type of the current cell
 */
Type extract_type(Options options);

/**
 * Function to swap old grid and the new one to continue the simulation with the same grid.
 * @param old old grid of the current timestep
 * @param new new grid of the next timestep
 * @param size grid's side length
 */
void swap_grids(Cell* old, Cell* new, int size);

#endif