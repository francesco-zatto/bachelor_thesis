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
int cells_B_number = 100;

/**
 * Default number of T helper lymphocytes of the simulation.
 */
int cells_T_number = 100;

/**
 * Default number of antigens of the simulation.
 */
int ag_number = 500;

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
 * @param positions positions taken by cells inside the grid
 * @param n number of cells in the simulation
 */
void save_positions(Vector* positions, int n);

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

    Vector* positions = calloc(options.total_number_cells, sizeof(Vector));
    for (int i = 0; i < options.total_number_cells; i++)
    {
        positions[i] = cells[i].position;
    }
    save_positions(positions, options.total_number_cells);

    free(positions);
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

void save_positions(Vector* positions, int n) 
{
    FILE* file = fopen("start.csv", "w");
    fprintf(file, "\n");
    for (int i = 0; i < n; i++) 
    {
        fprintf(file, "%f;%f;\n", positions[i].x, positions[i].y);
    }
    fprintf(file, "\n");
    fclose(file);
}