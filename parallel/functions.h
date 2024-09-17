#ifndef _FUNCTIONS_HEADER
#define _FUNCTIONS_HEADER

#include <stdbool.h>

#include "../cell.h"

/**
 * Fluid viscosity. It slows down the cell's movement in the fluid.
 */
#define LAMBDA 0.1

/**
 * Size of each misuration timestep.
 */
#define TIMESTEP 0.2

/**
 * Given a position, it returns a pointer to the requested cell of the host's grid.
 * @param grid grid of the simulation
 * @param position position (x, y) of the requested cell
 * @returns a pointer to a cell
 */
__host__ Cell* host_access_grid(Grid* grid, Vector position);

/**
 * Given a position, it returns a pointer to the requested cell of the device's grid.
 * @param grid grid of the simulation
 * @param position position (x, y) of the requested cell
 * @returns a pointer to a cell
 */
__device__ Cell* device_access_grid(Grid* grid, Vector position);

/**
 * It creates a duplicate in a certain position of a given B lymphocyte.
 * @param old_cell B cell to duplicate
 * @param duplicate duplicate cell to return as output parameter
 * @param position position to place the new duplicate
 */
__device__ void create_duplicate(Cell old_cell, Cell* duplicate, Vector position);

/**
 * It creates an antibody when a B cell is operative and ready to create them.
 * @param B_cell B lymphocyte that creates antigens
 * @param antibody antibody to return as output parameter
 * @param position position to place the new antibody
 */
__device__ void create_antibody(Cell B_cell, Cell* antigen, Vector position);


/**
 * It copies the value of a receptor in another. It can be used in cell duplication.
 * @param new_receptor char array as output parameter
 * @param old_receptor char array as input parameter
 */
__device__ void copy_receptor(unsigned char new_receptor[RECEPTOR_SIZE], unsigned char old_receptor[RECEPTOR_SIZE]);

/**
 * Action made by a B lymphocyte every timestep. Depending on its status, it checks the presence of
 * antigens or T helper lymphocytes nearby or it duplicates and create more antibodies.
 * @param b this lymphocyte that calls the action
 * @param old_grid grid where the B lymphocyte is moving
 * @param new_grid grid where the B lymphocyte will be in next timestep
 */
__device__ void lympho_B_action(Cell* b, Grid* old_grid, Grid* new_grid);

/**
 * Default action made by any other cell without a specific action like B lyphocytes or antibodies.
 * @param cell this cell that calls the action
 * @param old_grid grid where the cell is moving
 * @param new_grid grid where the cell will be in next timestep
 */
__device__ void default_action(Cell* cell, Grid* old_grid, Grid* new_grid);

/**
 * Action made by an inactive B lymphocyte or an antibody to search for antigens nearby.
 * @param cell this lymphocyte or this antibody that calls the action
 * @param old_grid grid where the cells are moving in this timestep
 * @param new_grid grif where the cells will move in next timestep
 */
__device__ void search_antigens(Cell* cell, Grid* old_grid, Grid* new_grid);

/**
 * Action made by an active B lymphocyte that is already bound to an antigen, so it searches for the closest lympho T
 * to react with, so that the B cell can become operative.
 * @param b B lymphocyte that searches in the grid
 * @param old_grid grid of cells where the B lymphocyte has to search
 */
__device__ void search_lympho_T(Cell* b, Grid* old_grid);

/**
 * It corrects a position that could be outside the grid.
 * @param position cell position to check
 * @param size of the grid
 */
__device__ void correct_position(Vector* position, int size);

/**
 * It returns true if the other cell is an antigen with a matching receptor.
 * @param cell antibody or B lymphocyte that is searching for antigens
 * @param other possible to antigen to check if it has a matching receptor
 */
__device__ bool is_matching_antigen(Cell cell, Cell other);

/**
 * It returns the hamming distance of two sequence of bits as array of chars, that are receptors.
 * @param receptor_cell first bit sequence
 * @param receptor_other second bit sequence
 */
__device__ int hamming_distance(unsigned char receptor_cell[RECEPTOR_SIZE], unsigned char receptor_other[RECEPTOR_SIZE]);

/**
 * Action made by an inactive B lymphocyte or an antibody when it finds an antigen nearby.
 * An antibody destroys it, freeing that cell, instead a B lymphocyte becomes active 
 * searching for a T lymphocyte.
 * @param cell this lymphocyte or this antibody that calls the action
 * @param other cell of the found antigen, so that the antibody can destroy it and free it
 */
__device__ void find_antigen(Cell* cell, Cell* other);

/**
 * It duplicates the given cell in a free nearby position in the new_grid.
 * The old grid and the new grid can be a reference to the same grid.
 * @param cell the cell to duplicate
 * @param old_grid the old grid where the cell was
 * @param new_grid the new grid where the cell's duplicate will be
 */
__device__ void duplicate(Cell* cell, Grid* old_grid, Grid* new_grid);

/**
 * Given an operative B lymphocyte, it creates antibodies with the same receptor nearby to destroy antigens.
 * @param cell the B lymphocyte that creates antibodies
 * @param old_grid the old grid where the cell was
 * @param new_grid the new grid where the created antibodies will be
 */
__device__ void create_antibodies(Cell* cell, Grid* old_grid, Grid* new_grid);

#endif