#ifndef _SIM_UTILS_HEADER
#define _SIM_UTILS_HEADER

#include "cell.h"

/**
 * Function to generate the grid, taking its sides.
 * @param grid grid to generate
 * @param length length of the square grid size
 * @param options generation options, like B cells quantity
 */
void generate(Cell* grid, int length, Options options);

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