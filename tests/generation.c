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

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>

#include "../cell.h"

/**
 * Size of grid size.
 */
#define SIZE 2000

/**
 * Default number of B lymphocytes of the simulation.
 */
int CELLS_B_NUMBER = 1000;

/**
 * Default number of T helper lymphocytes of the simulation.
 */
int CELLS_T_NUMBER = 1000;

/**
 * Default number of antigens of the simulation.
 */
int AG_NUMBER = 5000;

/**
 * Generation process of the simulation cells given the options.
 * @param cells array of cells to generate
 * @param options options with info about the generation of cells
 */
void generation(Cell* cells, Options options);

/**
 * It creates a cell of a certain given type.
 * @param cell cell pointer to create
 * @param type cell type
 */
void create_cell(Cell* cell, Type type);

/**
 * It returns true if the passed positions are equal, otherwise false.
 * @param a first position to compare
 * @param b second position to compare
 */
bool compare_positions(Vector a, Vector b);

/**
 * It save cells positions in an external file.
 * @param cells cells inside the grid
 * @param n number of cells in the simulation
 */
void save_positions(Cell* cells, int n);

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    if (argc == 4) 
    {
        CELLS_B_NUMBER = atoi(argv[1]);
        CELLS_T_NUMBER = atoi(argv[2]);
        AG_NUMBER = atoi(argv[3]);
    }

    Options options = 
    {
        .total_number_cells = CELLS_B_NUMBER + CELLS_T_NUMBER + AG_NUMBER,
        .cells_B_number = CELLS_B_NUMBER,
        .cells_T_number = CELLS_T_NUMBER,
        .ag_number = AG_NUMBER
    };
    Cell* cells = calloc(options.total_number_cells, sizeof(Cell));

    generation(cells, options);

    save_positions(cells, options.total_number_cells);

    free(cells);
    return 0;
}

void generation(Cell* cells, Options options)
{
    int cells_quantities[3] = 
    {
        options.cells_B_number, 
        options.cells_B_number + options.cells_T_number, 
        options.total_number_cells
    };
    int i = 0, actual = 0;
    Type current_type = B;
    while(i < options.total_number_cells) 
    {
        if (cells_quantities[current_type] == i)
        {
            current_type++;
        }
        create_cell(&cells[i], current_type);
        bool found = false;
        for (int j = 0; j < i; j++) 
        {
            if (compare_positions(cells[j].position, cells[i].position))
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            i++;
        }
    }
}

void create_cell(Cell* cell, Type type)
{
    cell->type = type;
    for (int i = 0; i < RECEPTOR_SIZE; i++) 
    {
        cell->receptor[i] = (char)(rand() % UCHAR_MAX);
    }
    cell->position.x = rand() % SIZE;
    cell->position.y = rand() % SIZE;     
}

bool inline compare_positions(Vector a, Vector b)
{
    return a.x == b.x && a.y == b.y;
}

void save_positions(Cell* cells, int n) 
{
    FILE* file = fopen("start.csv", "w");
    fprintf(file, "Type;X;Y\n");
    for (int i = 0; i < n; i++) 
    {
        fprintf(file, "%d;", cells[i].type);
        fprintf(file, "%f;%f\n", cells[i].position.x, cells[i].position.y);
    }
    fclose(file);
}