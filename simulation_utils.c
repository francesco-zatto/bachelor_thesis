#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "simulation_utils.h"
#include "functions.h"

int counts[4] = {0, 0, 0, 0};

void read_parameters(Options *options, int *parameters, int n)
{
    Options DEFAULT_OPTIONS = 
    {
    .cells_B_number = cells_B_number,
    .cells_T_number = cells_T_number,
    .grid_size = size,
    .ag_number = ag_number,
    .total_number_cells = cells_B_number + cells_T_number + ag_number
    };

    *options = DEFAULT_OPTIONS;
    switch (n)
    {
    case 4:
        options->grid_size = parameters[3];
    case 3:
        options->cells_B_number = parameters[0];
        options->cells_T_number = parameters[1];
        options->ag_number = parameters[2];
        break;
    case 1:
        options->grid_size = parameters[0];
        break;
    }

    int grid = options->grid_size * options->grid_size;
    int total_cells = options->ag_number + options->cells_B_number + options->cells_T_number;

    assert(grid < total_cells);
}

void generation(Cell *grid, int length, Options options)
{
    for (int i = 0; i < length; i++) 
    {
        for (int j = 0; j < length; j++)
        {
            Vector position = {i, j};
            Type type = extract_type(options);
            create_cell(access_grid(grid, position), position, type);
        }
    }
    printf("B: %d\tT: %d\t Ag: %d\t Free: %d\n", counts[B], counts[T], counts[Ag], counts[FREE]);
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
    double cells_B_prob = options.cells_B_number / total;
    double cells_T_prob = cells_B_prob + options.cells_T_number / total;
    double ag_prob = cells_T_prob + options.ag_number / total;

    double prob = rand() / RAND_MAX;
    Type type = prob < cells_B_prob 
            ? B 
            : (prob < cells_T_prob ? T : (prob < ag_prob ? Ag : FREE));
    counts[type]++;
    return type;
}
