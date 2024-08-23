#ifndef _FUNCTIONS_HEADER
#define _FUNCTIONS_HEADER

#include <stdbool.h>

#include "cell.h"

/**
 * Given a position, it returns a pointer to the requested cell of the grid.
 * @param grid grid of the simulation
 * @param position position (x, y) of the requested cell
 * @returns a pointer to a cell
 */
inline Cell* access_grid(Cell* grid, Vector position);

/**
 * Action made by a B lymphocyte every timestep. Depending on its status, it checks the presence of
 * antigens or T helper lymphocytes nearby or it duplicates and create more antibodies.
 * @param b this lymphocyte that calls the action
 * @param old_grid grid where the B lymphocyte is moving
 * @param new_grid grid where the B lymphocyte will be in next timestep
 */
void lympho_B_action(Cell* b, Cell* old_grid, Cell* new_grid);

/**
 * Action made by an inactive B lymphocyte or an antibody to search for antigens nearby.
 * @param cell this lymphocyte or this antibody that calls the action
 * @param grid grid where the cells are moving
 */
void search_antigens(Cell* cell, Cell* grid);

/**
 * It corrects a position that could be outside the grid.
 * @param position cell position to check
 */
void correct_position(Vector* position);

/**
 * It returns true if the other cell is an antigen with a matching receptor.
 * @param cell antibody or B lymphocyte that is searching for antigens
 * @param other possible to antigen to check if it has a matching receptor
 */
bool is_matching_antigen(Cell cell, Cell other);

/**
 * It returns the hamming distance of two sequence of bits as array of chars, that are receptors.
 * @param receptor_cell first bit sequence
 * @param receptor_other second bit sequence
 */
int hamming_distance(char receptor_cell[RECEPTOR_SIZE], char receptor_other[RECEPTOR_SIZE]);

/**
 * Action made by an inactive B lymphocyte or an antibody when it finds an antigen nearby.
 * An antibody destroys it, freeing that cell, instead a B lymphocyte becomes active 
 * searching for a T lymphocyte.
 * @param cell this lymphocyte or this antibody that calls the action
 * @param other cell of the found antigen, so that the antibody can destroy it and free it
 */
void find_antigen(Cell* cell, Cell* other);

/**
 * It duplicates the given cell in a free nearby position in the new_grid.
 * The old grid and the new grid can be a reference to the same grid.
 * @param cell the cell to duplicate
 * @param old_grid the old grid where the cell was
 * @param new_grid the new grid where the cell's duplicate will be
 */
void duplicate(Cell* cell, Cell* old_grid, Cell* new_grid);

/**
 * Given an operative B lymphocyte, it creates antibodies with the same receptor nearby to destroy antigens.
 * @param cell the B lymphocyte that creates antibodies
 * @param old_grid the old grid where the cell was
 * @param new_grid the new grid where the created antibodies will be
 */
void create_antibodies(Cell* cell, Cell* old_grid, Cell* new_grid);

#endif