#ifndef _FUNCTIONS_HEADER
#define _FUNCTIONS_HEADER

#include <stdbool.h>

#include "../cell.h"


/**
 * Action made by a B lymphocyte every timestep. Depending on its status, it checks the presence of
 * antigens or T helper lymphocytes nearby or it duplicates and create more antibodies.
 * @param b this lymphocyte that calls the action
 * @param grid grid where the B lymphocyte is moving
 */
void lympho_B_action(Cell* b, Cell* grid);

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

#endif