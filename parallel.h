#ifndef _PARALLEL_HEADER
#define _PARALLEL_HEADER

#include <stdbool.h>
#include <curand.h>
#include <curand_kernel.h>

#include "cuda_cell.h"

/**
 * Fluid viscosity. It slows down the cell's movement in the fluid.
 */
#define LAMBDA 0.1

/**
 * Size of each misuration timestep.
 */
#define TIMESTEP 0.2

#define PI 3.14159

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
 * @param size grid's size
 * @returns a pointer to a cell
 */
__device__ Cell* device_access_grid(Cell* grid, Vector position, int size);

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
 * @param size grid's size
 */
__device__ void lympho_B_action(Cell* b, Cell* old_grid, Cell* new_grid, int size);

/**
 * Default action made by any other cell without a specific action like B lyphocytes or antibodies.
 * @param cell this cell that calls the action
 * @param old_grid grid where the cell is moving
 * @param new_grid grid where the cell will be in next timestep
 * @param size grid's size
 */
__device__ void default_action(Cell* cell, Cell* old_grid, Cell* new_grid, int size);

/**
 * Action made by an inactive B lymphocyte or an antibody to search for antigens nearby.
 * @param cell this lymphocyte or this antibody that calls the action
 * @param old_grid grid where the cells are moving in this timestep
 * @param new_grid grif where the cells will move in next timestep
 * @param size grid's size
 */
__device__ void search_antigens(Cell* cell, Cell* old_grid, Cell* new_grid, int size);

/**
 * Action made by an active B lymphocyte that is already bound to an antigen, so it searches for the closest lympho T
 * to react with, so that the B cell can become operative.
 * @param b B lymphocyte that searches in the grid
 * @param old_grid grid of cells where the B lymphocyte has to search
 * @param size grid's size
 */
__device__ void search_lympho_T(Cell* b, Cell* old_grid, int size);

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
 * @param size grid's size
 */
__device__ void duplicate(Cell* cell, Cell* old_grid, Cell* new_grid, int size);

/**
 * Given an operative B lymphocyte, it creates antibodies with the same receptor nearby to destroy antigens.
 * @param cell the B lymphocyte that creates antibodies
 * @param old_grid the old grid where the cell was
 * @param new_grid the new grid where the created antibodies will be
 * @param size grid's size
 */
__device__ void create_antibodies(Cell* cell, Cell* old_grid, Cell* new_grid, int size);

/**
 * Function to compute acceleration using Langevin equation depending on given parameters.
 * @param velocity body current velocity
 * @param collision_forces collision forces felt by the body
 * @param mass body mass
 */
__device__ float langevin_equation(float velocity, float collision_forces, float mass);

/**
 * Getter for body mass depending on its type. For example, lymphocytes are way heavier than antibodies.
 * @param type cell type
 * @return cell's mass
 */
__device__ float get_mass(Type type);

/**
 * Using Langevin equation and following a brownian motion inside a fluid, it updates cells coordinates
 * inside 2D grid, so moving with discrete steps.
 * @param cell cell to update their current velocity and position
 * @param new_grid grid to place cell in next timestep
 * @param rand_state CUDA object that generates random numbers
 */
__device__ void movement(Cell *cell, Grid *new_grid, curandState* rand_state);

/**
 * It generates two random numbers following a normal distribution.
 * @param box_muller_number two numbers following a normal distribution generated by Box-Muller method.
 */
__device__ void box_muller(float box_muller_number[2], curandState* rand_state);

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
