#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <curand.h>
#include <curand_kernel.h>

#include "parallel/functions.h"
#include "parallel/physics.h"
#include "parallel/simulation_utils.h"

__host__ inline Cell* host_access_grid(Grid* grid, Vector position)
{
    /**
     * Because of the choice of using a single array grid, accessing it is way much more complicated,
     * but next line has the same meaning of &matrix[position.x][position.y] for a multi-array grid.
     */
    return &grid->matrix[(int)position.x * grid->size + (int)(position).y];
}

__device__ inline Cell* device_access_grid(Grid* grid, Vector position)
{
    /**
     * Same as host_access_grid function, but in device memory
     */
    return &grid->matrix[(int)position.x * grid->size + (int)(position).y];
}

__device__ void lympho_B_action(Cell* b, Grid* old_grid, Grid* new_grid)
{
    switch (b->status)
    {
        //If inactive, it searches for antigens.
        case INACTIVE: 
        {
            search_antigens(b, old_grid, new_grid);
            break;
        }
        //If active, it searches for a T lymphocytes.
        case ACTIVE:
        {
            search_lympho_T(b, old_grid);
        }
        //If operative, it duplicates itself and create antibodies.
        case OPERATIVE:
        {
            duplicate(b, old_grid, new_grid);
            create_antibodies(b, old_grid, new_grid);
            break;
        }
    }
}

__device__ void default_action(Cell *cell, Grid *old_grid, Grid *new_grid)
{
    //Empty action for entities with any particular action to do.
}

__device__ void search_antigens(Cell* cell, Grid* old_grid, Grid* new_grid)
{
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE; j++) 
        {
            /**
             * Checking if the position has an antigen with a matching antigen.
             * In that case, the B cell has found it and it change its status.
             */
            Vector current_position = 
            {
                .x = cell->position.x + i,
                .y = cell->position.y + j 
            };
            correct_position(&current_position, old_grid->size);
            Cell* other = NULL; //device_access_grid(old_grid, current_position);
            if (is_matching_antigen(*cell, *other))
            {
                find_antigen(cell, other);
                return;
            }
        }
    }
}

__device__ bool is_matching_antigen(Cell cell, Cell other) 
{
    //Hamming distance of the two receptors has to be greater or equal than the threshold.
    return other.type == Ag && hamming_distance(cell.receptor, other.receptor) >= AFFINITY_MIN;
}

__device__ int hamming_distance(unsigned char receptor_cell[RECEPTOR_SIZE], unsigned char receptor_other[RECEPTOR_SIZE])
{
    int distance = 0;
    char xor_receptor[RECEPTOR_SIZE] = {0};
    /**
     * Computing xor between every char of the two receptors' arrays.
     */
    for (int i = 0; i < RECEPTOR_SIZE; i++)
    {
        xor_receptor[i] = (receptor_cell[i] ^ receptor_other[i]);
    }
    /**
     * Computing the hamming distance as the number of the bit set to 1.
     */
    for (int i = 0; i < RECEPTOR_SIZE; i++) 
    {
        for (int j = 0; j < 8; j++)
        {
            distance += (xor_receptor[i] >> j) & 0x1;
        }
    }
    return distance;
}

__device__ void correct_position(Vector* position, int size) 
{
    if (position->x >= size)
    {
        position->x -= size;
    }
    if (position->x < 0)
    {
        position->x = size + position->x;
    }
    if (position->y >= size)
    {
        position->y -= size;
    }
    if (position->y < 0)
    {
        position->y = size + position->y;
    }
}

__device__ void find_antigen(Cell* cell, Cell* other)
{
    switch (cell->type)
    {
        //A B cell becomes active when it finds an antigen.
        case B:
        {
            cell->status = ACTIVE;
            break;
        }
        //An antibody kills every antigen that it finds.
        case Ab:
        {
            other->type = FREE;
            break;
        }
    }
}

__device__ void search_lympho_T(Cell* b, Grid* old_grid)
{
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE; j++) 
        {
            /**
             * Checking if the position is taken by a T lymphocyte.
             * In that case, the B cell becomes operative.
             */
            Vector current_position = 
            {
                .x = b->position.x + i,
                .y = b->position.y + j 
            };
            correct_position(&current_position, old_grid->size);
            Cell* other = NULL; //device_access_grid(old_grid, current_position);
            if (other->type == T)
            {
                b->status = OPERATIVE;
                return;
            }
        }
    }
}

__device__ void duplicate(Cell* cell, Grid* old_grid, Grid* new_grid)
{
    bool duplicated = false;
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE && !duplicated; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE && !duplicated; j++) 
        {
            /**
             * When the B cell finds a free cell nearby it duplicates.
             * In that case it creates a duplicate and becomes inactive.
             */
            Vector new_position = 
            {
                .x = cell->position.x + i,
                .y = cell->position.y + j
            };
            correct_position(&new_position, old_grid->size);
            Cell* free_cell = NULL; //device_access_grid(new_grid, new_position);
            if (free_cell->type == FREE)
            {
                create_duplicate(*cell, free_cell, new_position);
                cell->status = INACTIVE;
                duplicated = true;
            }
        }
    }
}

__device__ void create_duplicate(Cell old_cell, Cell* duplicate, Vector position)
{
    duplicate->action = old_cell.action;
    duplicate->type = old_cell.type;
    duplicate->status = INACTIVE;
    duplicate->velocity.x = 0;
    duplicate->velocity.y = 0;
    duplicate->position = position;
    copy_receptor(duplicate->receptor, old_cell.receptor);
}

__device__ void copy_receptor(unsigned char new_receptor[RECEPTOR_SIZE], unsigned char old_receptor[RECEPTOR_SIZE])
{
    for (int i = 0; i < RECEPTOR_SIZE; i++)
    {
        new_receptor[i] = old_receptor[i];
    }
}

__device__ void create_antibodies(Cell* cell, Grid* old_grid, Grid* new_grid)
{
    int created = 0;
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE && created < NUMBER_CREATED_ANTIGENS; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE && created < NUMBER_CREATED_ANTIGENS; j++) 
        {
            /**
             * Checking it the position is free, then it creates an antibody in that position.
             */
            Vector new_position = 
            {
                .x = cell->position.x + i,
                .y = cell->position.y + j
            };
            correct_position(&new_position, old_grid->size);
            Cell* free_cell = NULL; device_access_grid(new_grid, new_position);
            if (free_cell->type == FREE)
            {
                create_antibody(*cell, free_cell, new_position);
                created++;
            }
        }
    }
}

__device__ void create_antibody(Cell B_cell, Cell* antibody, Vector position)
{
    antibody->action = search_antigens;
    antibody->type = Ab;
    antibody->velocity.x = 0;
    antibody->velocity.y = 0;
    antibody->position = position;
    copy_receptor(antibody->receptor, B_cell.receptor);
}


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
                //TODO KEEP correct_position(&position, grid->size);
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

/**
 * Individual timestep ran by each thread, it will be the main kernel of the simulation.
 * @param grid grid of the current iteration
 * @param next_grid grid of the next iteration
 */
__global__ void timestep(Grid* grid, Grid* next_grid, curandState* rand_state)
{
    /**
     * Selecting a cell in position (x, y), executing its action and making it move.
     */
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < grid->size  && y < grid->size)
    {
        Vector position = {(float) x,(float) y};
        Cell* cell = NULL; //device_access_grid(grid, position);
        cell->action(cell, grid, next_grid);  
        //movement(cell, next_grid, rand_state);
    }
}

/**
 * Actual simulation function that for n timesteps updates cells in the grid.
 * @param grid grid of the current iteration
 * @param next_grid grid of the next iteration
 * @param options simulation options
 */
__host__ void simulation(Grid* grid, Grid* next_grid, Options options)
{
    dim3 block_size = {32, 32, 1};
    unsigned int grid_length_x = (grid->size) / block_size.x + 1;
    unsigned int grid_length_y = (grid->size) / block_size.y + 1;
    dim3 number_blocks = {grid_length_x, grid_length_y, 1};
    dim3 thread_per_block = {block_size.x, block_size.y, 1};

    curandState *rand_state;
    cudaMalloc(&rand_state, sizeof(curandState));

    for (int t = 0; t < TIMESTEPS; t++)
    {
        timestep<<<number_blocks, thread_per_block>>>(grid, next_grid, rand_state);
        /**
         * Swapping grids, so that for every iteration it is used the same grid variable.
         */
        swap_grids<<<number_blocks, thread_per_block>>>(grid, next_grid, options.grid_size);
    }

    cudaFree(rand_state);
}

int main(int argc, char const *argv[])
{
    srand(time(NULL));

    /**
     * Initialization of simulation options and allocation of grids' memory.
     */
    Options options;
    read_parameters(&options, argv, argc - 1);

    Grid h_grid, h_next_grid,
        d_grid, d_next_grid;

    int malloc_size = sizeof(Cell) * options.grid_size * options.grid_size;

    d_grid.size = d_next_grid.size = h_grid.size = h_next_grid.size = options.grid_size;
    h_grid.matrix = (Cell*)malloc(malloc_size),
    h_next_grid.matrix = (Cell*)malloc(malloc_size);
    cudaMalloc((void**) &(d_grid.matrix), malloc_size);
    cudaMalloc((void**) &(d_next_grid.matrix), malloc_size);

    /**
     * Placing the cells and antigens in the grid and placing only free cells in the next iteration grid.
     */
    generation(&h_grid, options);
    free_grid(&h_next_grid);
    cudaMemcpy((void*) d_grid.matrix, (void*) h_grid.matrix, malloc_size, cudaMemcpyHostToDevice);
    cudaMemcpy((void*) d_next_grid.matrix, (void*) h_next_grid.matrix, malloc_size, cudaMemcpyHostToDevice);

    //Save cells at the start in a file
    const char* start_file = "grids/start.csv"; 
    save_grid(&h_grid, start_file);

    simulation(&d_grid, &d_next_grid, options);

    //Save cells in the middle of the simulation, before inserting new antigens.
    cudaMemcpy((void*) h_grid.matrix, (void*) d_grid.matrix, malloc_size, cudaMemcpyDeviceToHost);
    const char* mid_file = "grids/mid.csv";
    save_grid(&h_grid, mid_file);
    insert_antigens(&h_grid);
    cudaMemcpy((void*) d_grid.matrix, (void*) h_grid.matrix, malloc_size, cudaMemcpyHostToDevice);

    simulation(&d_grid, &d_next_grid, options);

    //Save cells at the end in a file
    cudaMemcpy((void*) h_grid.matrix, (void*) d_grid.matrix, malloc_size, cudaMemcpyDeviceToHost);
    const char* end_file = "grids/end.csv";
    save_grid(&h_grid, end_file);

    free(h_grid.matrix);
    free(h_next_grid.matrix);
    cudaFree(d_grid.matrix);
    cudaFree(d_next_grid.matrix);
    return 0;
}