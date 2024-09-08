#ifndef _SIM_UTILS_HEADER
#define _SIM_UTILS_HEADER

#include "cell.h"

/**
 * Size of grid size.
 */
const int size = 2000;

/**
 * Default number of B lymphocytes of the simulation.
 */
const int cells_B_number = 1000;

/**
 * Default number of T helper lymphocytes of the simulation.
 */
const int cells_T_number = 1000;

/**
 * Default number of antigens of the simulation.
 */
const int ag_number = 5000;

/**
 * Function to read command line parameters to change default values,
 * like numbers of cells or grid size.
 * @param options output with the setted options
 * @param parameters command line parameters
 * @param n number of parameters to save
 */
void read_parameters(Options* options, int* parameters, int n);

/**
 * Function to generate the grid, taking its sides.
 * @param grid grid to generate
 * @param length length of the square grid size
 * @param options generation options, like B cells quantity
 */
void generation(Cell* grid, int length, Options options);

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

#endif