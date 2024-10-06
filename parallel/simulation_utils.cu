#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "simulation_utils.h"
#include "functions.h"

const Vector NULL_VECTOR = {0, 0};

__constant__ Cell D_FREE_CELL = 
{
    .type = FREE,
    .action = default_action
};

const Cell H_FREE_CELL = 
{
    .type = FREE,
    .action = default_action
};

const int TIMESTEPS = 100;

const int GRID_SIZE = 500;

const int CELLS_B_NUMBER = 200;

const int CELLS_T_NUMBER = 200;

const int AG_NUMBER = 2000;

const int NEW_AG_NUMBER = 8000;

__host__ void read_parameters(Options *options, const char *parameters[], int n)
{
    Options DEFAULT_OPTIONS = 
    {
        .total_number_cells = CELLS_B_NUMBER + CELLS_T_NUMBER + AG_NUMBER,
        .cells_B_number = CELLS_B_NUMBER,
        .cells_T_number = CELLS_T_NUMBER,
        .ag_number = AG_NUMBER,
        .grid_size = GRID_SIZE
    };

    //Setting simulation options to default.
    *options = DEFAULT_OPTIONS;
    switch (n)
    {
        //For 4 user parameters, the fourth is the grid's size.
        case 4:
            options->grid_size = atoi(parameters[4]);
        //First 3 user parameters are B cells, T cells and antigens number.
        case 3:
            options->cells_B_number = atoi(parameters[1]);
            options->cells_T_number = atoi(parameters[2]);
            options->ag_number = atoi(parameters[3]);
            options->total_number_cells = options->cells_B_number + options->cells_T_number + options->ag_number;
            break;
        //In case of only 1 parameter, it's the grid's size.
        case 1:
            options->grid_size = atoi(parameters[1]);
            break;
    }
    assert(options->grid_size > sqrt(options->total_number_cells));
}

__host__ void generation(Grid *grid, Options options)
{
    for (int i = 0; i < grid->size; i++) 
    {
        for (int j = 0; j < grid->size; j++)
        {
            /**
             * Extracting type of the current position and creating a cell in it.
             */
            Vector position = {(float) i, (float) j};
            Type type = extract_type(options);
            create_cell(host_access_grid(grid, position), position, type);
        }
    }
}

__host__ void create_cell(Cell *cell, Vector position, Type type)
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

__host__ Type extract_type(Options options)
{
    /**
     * It chooses a type as a number in [0, size * size).
     * For example, if n is in [B_cells, B_cells + T_cells] it extract a T cell.
     * So, it doesn't extract exactly the requested number of cells, but a very close number.
     */
    int total = options.grid_size * options.grid_size;
    int cells_B_prob = options.cells_B_number;
    int cells_T_prob = cells_B_prob + options.cells_T_number;
    int ag_prob = cells_T_prob + options.ag_number;

    int prob = rand() % total;
    Type type = prob < cells_B_prob 
            ? B 
            : (prob < cells_T_prob ? T : (prob < ag_prob ? Ag : FREE));
    return type;
}

__global__ void swap_grids(Grid *old_grid, Grid *new_grid, int size)
{
    Cell *old_cell, *new_cell;
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < size && y < size)
    {
        Vector position = {(float) x, (float) y};
        old_cell = NULL; //device_access_grid(old_grid, position);
        new_cell = NULL; //device_access_grid(new_grid, position);
        *old_cell = *new_cell;
        *new_cell = D_FREE_CELL;
    }
}

__host__ void free_grid(Grid* grid)
{
    for (int i = 0; i < grid->size; i++)
    {
        for (int j = 0; j < grid->size; j++)
        {
            Vector position = {(float) i, (float) j};
            Cell* cell = host_access_grid(grid, position);
            *cell = H_FREE_CELL;
        }
    }
}

__host__ void save_grid(Grid* grid, const char* filename)
{
    FILE* out = fopen(filename, "w");

    //Header of the csv file
    fprintf(out, "Type;Pos_x;Pos_y;Vel_x;Vel_y;");
    for (int l = 0; l < RECEPTOR_SIZE; l++)
    {
        fprintf(out, "Receptor_%d;", l);
    }
    fprintf(out, "Status\n");

    //A row for each cell and a column for each cell's member.
    for (int i = 0; i < grid->size; i++)
    {
        for (int j = 0; j < grid->size; j++)
        {
            Vector position = {(float) i, (float) j};
            Cell cell = *host_access_grid(grid, position);
            if (cell.type != FREE)
            {
                fprintf
                (
                        out, "%d;%f;%f;%f;%f;", 
                        cell.type, cell.position.x, cell.position.y, cell.velocity.x, cell.velocity.y
                );
                for (int l = 0; l < RECEPTOR_SIZE; l++)
                {
                    fprintf(out, "%d;", cell.receptor[l]);
                }
                fprintf(out, "%d\n", cell.status);
            }
        }
    }
    fclose(out);
}

__host__ void insert_antigens(Grid *grid)
{
    int step = grid->size / 10;
    int inserted = 0;
    for (int start = 0; start < step; start++)
    {
        for (int i = start; i < grid->size; i += step)
        {
            for (int j = start; j < grid->size; j += step)
            {
                /**
                 * Checking if current position is free, in that case insert an antigen.
                 * If the number of new antigen is already reached, then return.
                 */
                Vector position = {(float) i, (float) j};
                correct_position(&position, grid->size);
                Cell* cell = host_access_grid(grid, position);
                if (cell->type == FREE)
                {
                    create_cell(cell, position, Ag);
                    inserted++;
                }
                if (inserted == NEW_AG_NUMBER)
                {
                    return;
                }
            }
        }
    }
}
