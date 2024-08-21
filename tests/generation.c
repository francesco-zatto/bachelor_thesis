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
int cells_B_number = 1000;

/**
 * Default number of T helper lymphocytes of the simulation.
 */
int cells_T_number = 1000;

/**
 * Default number of antigens of the simulation.
 */
int ag_number = 5000;

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
        cells_B_number = atoi(argv[1]);
        cells_T_number = atoi(argv[2]);
        ag_number = atoi(argv[3]);
    }

    Options options = 
    {
        .total_number_cells = cells_B_number + cells_T_number + ag_number,
        .cells_B_number = cells_B_number,
        .cells_T_number = cells_T_number,
        .ag_number = ag_number
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