#include <stdlib.h>
#include <stdbool.h>

#include "simulation_utils.h"
#include "functions.h"

int counts[4] = {0, 0, 0, 0};

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
