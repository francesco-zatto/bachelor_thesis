#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "simulation_utils.h"
#include "functions.h"

int counts[5] = {0};

const Vector NULL_VECTOR = {0, 0};

const Cell FREE_CELL = {.type = FREE};

const int TIMESTEPS = 1000;

const int GRID_SIZE = 2000;

const int CELLS_B_NUMBER = 1000;

const int CELLS_T_NUMBER = 1000;

const int AG_NUMBER = 5000;

void read_parameters(Options *options, const char *parameters[], int n)
{
    Options DEFAULT_OPTIONS = 
    {
    .cells_B_number = CELLS_B_NUMBER,
    .cells_T_number = CELLS_T_NUMBER,
    .grid_size = GRID_SIZE,
    .ag_number = AG_NUMBER,
    .total_number_cells = CELLS_B_NUMBER + CELLS_T_NUMBER + AG_NUMBER
    };

    *options = DEFAULT_OPTIONS;
    switch (n)
    {
    case 4:
        options->grid_size = atoi(parameters[4]);
    case 3:
        options->cells_B_number = atoi(parameters[1]);
        options->cells_T_number = atoi(parameters[2]);
        options->ag_number = atoi(parameters[3]);
        options->total_number_cells = options->cells_B_number + options->cells_T_number + options->ag_number;
        break;
    case 1:
        options->grid_size = atoi(parameters[1]);
        break;
    }
    assert(options->grid_size > sqrt(options->total_number_cells));
}

void generation(Cell *grid, Options options)
{
    for (int i = 0; i < options.grid_size; i++) 
    {
        for (int j = 0; j < options.grid_size; j++)
        {
            Vector position = {i, j};
            Type type = extract_type(options);
            counts[type]++;
            create_cell(access_grid(grid, position), position, type);
        }
    }
    printf("B: %d, T: %d, Ag: %d, Free: %d\n", counts[B], counts[T], counts[Ag], counts[FREE]);
}

void create_cell(Cell *cell, Vector position, Type type)
{
    cell->position = position;
    cell->velocity = NULL_VECTOR;
    cell->type = type;
    cell->status = INACTIVE;
    for (int i = 0; i < RECEPTOR_SIZE; i++) 
    {
        cell->receptor[i] = (char)(rand() % UCHAR_MAX);
    }
    cell->action = cell->type == B 
            ? lympho_B_action 
            : (cell->type == Ab ? search_antigens : default_action);
}

Type extract_type(Options options)
{
    int total = options.grid_size * options.grid_size;
    int cells_B_prob = options.cells_B_number;
    int cells_T_prob = cells_B_prob + options.cells_T_number;
    int ag_prob = cells_T_prob + options.ag_number;

    int prob = rand() % total;
    Type type = prob < cells_B_prob 
            ? B 
            : (prob < cells_T_prob ? T : (prob < ag_prob ? Ag : FREE));
    counts[type]++;
    return type;
}

void swap_grids(Cell *old, Cell *new, int size)
{
    Cell *old_cell, *new_cell;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            Vector position = {i, j};
            old_cell = access_grid(old, position);
            new_cell = access_grid(new, position);
            *old_cell = *new_cell;
            *new_cell = FREE_CELL;
        }
    }
}
