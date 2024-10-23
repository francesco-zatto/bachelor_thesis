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
void read_parameters(Options* options, const char* parameters[], int n);

/**
 * Function to generate the grid, taking generation's options.
 * @param grid grid to generate
 * @param options generation options, like B cells quantity
 */
void generation(Grid* grid, Options options);

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
 * @param old_grid old grid of the current timestep
 * @param new_grid new grid of the next timestep
 * @param size grid's side length
 */
void swap_grids(Grid* old_grid, Grid* new_grid, int size);

/**
 * Function to fill a grid with only free cells. It does not deallocate grid memory.
 * @param grid output grid that will contain only free cells
 */
void free_grid(Grid* grid);

/**
 * Function to save every information about the cells inside the given grid in the current state
 * @param grid grid to save in external file
 * @param filename name of the file
 */
void save_grid(Grid* grid, char* filename);

/**
 * Function to insert new antigens in the simulation's grid to test HIS memory.
 * @param grid grid which antigens will be inserted in
 */
void insert_antigens(Grid* grid);

#endif