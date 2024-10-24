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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include <curand.h>
#include <curand_kernel.h>

#include "parallel.h"
#include "cuda_cell.h"

#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

__host__
void host_correct_position(Vector* position, int size) 
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


__device__
void device_correct_position(Vector* position, int size) 
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

__host__
inline Cell* host_access_grid(Grid* grid, Vector position)
{
    /**
     * Because of the choice of using a single array grid, accessing it is way much more complicated,
     * but next line has the same meaning of &matrix[position.x][position.y] for a multi-array grid.
     */
    return &grid->matrix[(int)position.x * grid->size + (int)(position).y];
}

__device__
inline Cell* device_access_grid(Cell* grid, Vector position, int size)
{
    /**
     * Same as host_access_grid function, but in device memory
     */
    return &grid[(int)position.x * size + (int)(position).y];
}

/**
 * Function to find a free cell nearby if given starting cell is taken.
 * @param original original cell position in the grid, before the movement
 * @param start start cell to search for another free cell
 * @param grid grid where the search is taken
 * @param rand_state random number generator in CUDA
 * @param suze grid's size
 * @param return free cell near start
 */
__device__
static Cell* find_free_cell_nearby(Cell* original, Cell* start, Cell* grid, curandState* rand_state, int size)
{
    Cell* cell;
    int sign_x = (int)(start->position.x) % 2 == 0 ? +1: -1;
    int sign_y = (int)(start->position.y) % 2 == 0 ? +1: -1;
    int prox_dist = PROXIMITY_DISTANCE;
    for (int n = 0; n < prox_dist * prox_dist; n++)
    {
        /**
         * Looking for a nearby position and checking if its free. In that case, that cell will be taken.
         */
	float rand_n = curand_uniform(&rand_state[(int)(start->position.x * size + start->position.y)]);
	rand_n *= prox_dist;
	int i = (int)rand_n % prox_dist;
	int j = (int)(rand_n / prox_dist);
        Vector position = {start->position.x + sign_x * i, start->position.y + sign_y * j};
        device_correct_position(&position, size);
        cell = device_access_grid(grid, position, size);
        int old_type = atomicCAS((int*) &(cell->type), FREE, (int) original->type);
	if (old_type == FREE)
	{
            return cell;
        }
    }
    cell = device_access_grid(grid, original->position, size);
    atomicExch((int*) &(cell->type), (int) original->type);
    return cell;
}

__global__
void movement(Cell *old_grid, Cell *new_grid, curandState* rand_state, int size)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < size && y < size)
    {
        Vector position = {(float) x, (float) y};
        Cell* cell = device_access_grid(old_grid, position, size);
        //If cell is free, ignore it.
        if (cell->type != FREE)
	{
            /**
             * Computing box muller numbers for forces felt by the body and getting body mass.
             */
            float box_muller_number[2];
            box_muller_number[0] = curand_uniform(&(rand_state[x * size + y]));
	    box_muller_number[1] = curand_uniform(&(rand_state[x * size + y]));
	    float mass = get_mass(cell->type);

            /**
             * Computing deltaV with Langevin equation and then updating cell's position and velocity.
             */
            Vector delta_velocity =
            {
                .x = langevin_equation(cell->velocity.x, box_muller_number[0], mass),
                .y = langevin_equation(cell->velocity.y, box_muller_number[1], mass)
            };
	    Vector velocity =
	    {
	        .x = cell->velocity.x + delta_velocity.x,
	        .y = cell->velocity.y + delta_velocity.y
	    };
            Vector new_position =
	    {
	        .x = position.x + (float) round(velocity.x * TIMESTEP),
	        .y = position.y + (float) round(velocity.y * TIMESTEP)
	    };

           /**
            * Checking if the computed position is inside the grid and if it is free.
            */
            device_correct_position(&new_position, size);
            Cell* new_cell = device_access_grid(new_grid, new_position, size);
            int old_type = atomicCAS((int*) &(new_cell->type), FREE, cell->type);
            if (old_type != FREE)
            {
                new_cell = find_free_cell_nearby(cell, new_cell, new_grid, rand_state, size);
            }
            new_cell->position.x = new_position.x;
	    new_cell->velocity.x = velocity.x;
	    new_cell->position.y = new_position.y;
	    new_cell->velocity.y = velocity.y;
	    new_cell->status = cell->status;
	    for (int i = 0; i < RECEPTOR_SIZE; i++)
	    {
	    	new_cell->receptor[i] = cell->receptor[i];
	    }
	}
    }
}

__device__
float inline langevin_equation(float velocity, float collision_forces, float mass)
{
    return (-LAMBDA * velocity + collision_forces) / mass * TIMESTEP;
}

__device__
float inline get_mass(Type type)
{
    //Lymphocytes cells have a much higher mass than antigens and antibodies.
    switch (type)
    {
    case B:
    case T:
        return 0.2;
    case Ag:
    case Ab:
        return 0.01;
    default:
        assert(type == B || type == T || type == Ag || type == Ab);
    }
}

__device__
void box_muller(float box_muller_number[2], curandState* rand_state)
{
    float 
        u1 = (float)(curand(rand_state)) / (float)(RAND_MAX),
        u2 = (float)(curand(rand_state)) / (float)(RAND_MAX);
    box_muller_number[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
    box_muller_number[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2);
}

__device__
void lympho_B_action(Cell* b, Cell* old_grid, Cell* new_grid, int size)
{
    switch (b->status)
    {
        //If inactive, it searches for antigens.
        case INACTIVE: 
        {
            search_antigens(b, old_grid, new_grid, size);
            break;
        }
        //If active, it searches for a T lymphocytes.
        case ACTIVE:
        {
            search_lympho_T(b, old_grid, size);
	    break;
        }
        //If operative, it duplicates itself and create antibodies.
        case OPERATIVE:
        {
            duplicate(b, old_grid, new_grid, size);
            create_antibodies(b, old_grid, new_grid, size);
            break;
        }
    }
}

__device__
void default_action(Cell *cell, Cell *old_grid, Cell *new_grid, int size)
{
    //Empty action for entities with any particular action to do.
}

__device__
void search_antigens(Cell* cell, Cell* old_grid, Cell* new_grid, int size)
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
            device_correct_position(&current_position, size);
            Cell* other = device_access_grid(old_grid, current_position, size);
            if (is_matching_antigen(*cell, *other))
            {
                find_antigen(cell, other);
                return;
            }
        }
    }
}

__device__
bool is_matching_antigen(Cell cell, Cell other) 
{
    //Hamming distance of the two receptors has to be greater or equal than the threshold.
    return other.type == Ag && hamming_distance(cell.receptor, other.receptor) >= AFFINITY_MIN;
}

__device__
int hamming_distance(unsigned char receptor_cell[RECEPTOR_SIZE], unsigned char receptor_other[RECEPTOR_SIZE])
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

__device__
void find_antigen(Cell* cell, Cell* other)
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

__device__
void search_lympho_T(Cell* b, Cell* old_grid, int size)
{
    for (int i = +PROXIMITY_DISTANCE; i >= -PROXIMITY_DISTANCE; i--)
    {
        for (int j = +PROXIMITY_DISTANCE; j >= -PROXIMITY_DISTANCE; j--) 
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
            device_correct_position(&current_position, size);
            Cell* other = device_access_grid(old_grid, current_position, size);
            if (other->type == T)
            {
                b->status = OPERATIVE;
                return;
            }
        }
    }
}

__device__
void duplicate(Cell* cell, Cell* old_grid, Cell* new_grid, int size)
{
    bool duplicated = false;
    int sign_x = (int)(cell->position.x) % 10 >= 5 ? +1: -1;
    int sign_y = (int)(cell->position.y) % 10 >= 5 ? +1: -1;
    for (int i = sign_x * PROXIMITY_DISTANCE; i >= -sign_x * PROXIMITY_DISTANCE && !duplicated; i = sign_x > 0 ? i - 1 : i + 1)
    {
        for (int j = sign_y * PROXIMITY_DISTANCE; j >= -sign_y * PROXIMITY_DISTANCE && !duplicated; j = sign_y > 0 ? j - 1 : j + 1) 
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
            device_correct_position(&new_position, size);
            Cell* free_cell = device_access_grid(new_grid, new_position, size);
            int old_type = atomicCAS((int*) &(free_cell->type), FREE, cell->type);
	    if (old_type == FREE)
            {
                create_duplicate(*cell, free_cell, new_position);
                cell->status = INACTIVE;
                duplicated = true;
            }
        }
    }
}

__device__
void create_duplicate(Cell old_cell, Cell* duplicate, Vector position)
{
    duplicate->action = old_cell.action; 
    duplicate->type = old_cell.type;
    duplicate->status = INACTIVE;
    duplicate->velocity.x = 0;
    duplicate->velocity.y = 0;
    duplicate->position = position;
    copy_receptor(duplicate->receptor, old_cell.receptor);
}

__device__
void copy_receptor(unsigned char new_receptor[RECEPTOR_SIZE], unsigned char old_receptor[RECEPTOR_SIZE])
{
    for (int i = 0; i < RECEPTOR_SIZE; i++)
    {
        new_receptor[i] = old_receptor[i];
    }
}

__device__
void create_antibodies(Cell* cell, Cell* old_grid, Cell* new_grid, int size)
{
    int created = 0;
    int sign_x = (int)(cell->position.x) % 10 >= 5 ? +1 : -1;
    int sign_y = (int)(cell->position.y) % 10 >= 5 ? +1 : -1;
    for (int i = sign_x * PROXIMITY_DISTANCE; i >= -sign_x * PROXIMITY_DISTANCE && created < NUMBER_CREATED_ANTIGENS; i = sign_x > 0 ? i - 1 : i + 1)
    {
        for (int j = sign_y * PROXIMITY_DISTANCE; j >= -sign_y * PROXIMITY_DISTANCE && created < NUMBER_CREATED_ANTIGENS; j = sign_y > 0 ? j - 1 : j + 1) 
        {
            /**
             * Checking it the position is free, then it creates an antibody in that position.
             */
            Vector new_position = 
            {
                .x = cell->position.x + i,
                .y = cell->position.y + j
            };
            device_correct_position(&new_position, size);
            Cell* free_cell = device_access_grid(new_grid, new_position, size);
            int old_type = atomicCAS((int*) &(free_cell->type), FREE, Ab);
	    if (old_type == FREE)
            {
                create_antibody(*cell, free_cell, new_position);
                created++;
            }
        }
    }
}

__device__
void create_antibody(Cell B_cell, Cell* antibody, Vector position)
{
    antibody->action = search_antigens;
    antibody->type = Ab;
    antibody->velocity.x = 0;
    antibody->velocity.y = 0;
    antibody->position = position;
    copy_receptor(antibody->receptor, B_cell.receptor);
}


const Vector NULL_VECTOR = {0, 0};

__constant__
Cell D_FREE_CELL = 
{
    .type = FREE,
    .action = default_action
};

const Cell H_FREE_CELL = 
{
    .type = FREE,
    .position = NULL_VECTOR,
    .velocity = NULL_VECTOR,
    .status = INACTIVE,
    .action = default_action
};

const int TIMESTEPS = 100;

const int GRID_SIZE = 500;

const int CELLS_B_NUMBER = 200;

const int CELLS_T_NUMBER = 200;

const int AG_NUMBER = 2000;

const int NEW_AG_NUMBER = 500;

__host__
void read_parameters(Options *options, const char *parameters[], int n)
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

__host__
void generation(Grid *grid, Options options)
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

__host__
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

__host__
Type extract_type(Options options)
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

__global__
void swap_grids(Cell *old_grid, Cell *new_grid, int sub_grid_size, int grid_size, Vector start_position)
{
    Cell *old_cell, *new_cell;
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < sub_grid_size && y < sub_grid_size)
    {
        Vector position = {start_position.x + (float) x, start_position.y + (float) y};
        old_cell = device_access_grid(old_grid, position, grid_size);
        new_cell = device_access_grid(new_grid, position, grid_size);
	if (old_cell->type != FREE || new_cell->type != FREE)
	{
	    old_cell->type = new_cell->type;
	    old_cell->position.x = new_cell->position.x;
	    old_cell->velocity.x = new_cell->velocity.x;
	    old_cell->position.y = new_cell->position.y;
	    old_cell->velocity.y = new_cell->velocity.y;
	    old_cell->status = new_cell->status;
	    for (int i = 0; i < RECEPTOR_SIZE; i++)
	    {
		old_cell->receptor[i] = new_cell->receptor[i];
	    }
            new_cell->type = FREE;
	    new_cell->position.x = 0;
	    new_cell->velocity.x = 0;
	    new_cell->position.y = 0;
	    new_cell->velocity.y = 0;
	    new_cell->status = INACTIVE;
	    for (int i = 0; i < RECEPTOR_SIZE; i++)
	    {
		new_cell->receptor[i] = 0;
	    }
	}
    }
}

__host__
void clear_grid(Grid* grid)
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

__host__
void save_grid(Grid* grid, const char* filename)
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
                        cell.type, position.x, position.y, cell.velocity.x, cell.velocity.y
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

__host__
void insert_antigens(Grid *grid)
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
                host_correct_position(&position, grid->size);
                Cell* cell = host_access_grid(grid, position);
                if (cell->type == FREE)
                {
                    create_cell(cell, position, Ag);
                    inserted++;
                }
                if (inserted >= NEW_AG_NUMBER)
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
 * @param size grid's size
 */
__global__
void timestep(Cell* grid, Cell* next_grid, curandState* rand_state, int size)
{
    /**
     * Selecting a cell in position (x, y), executing its action and making it move.
     */
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < size  && y < size)
    {
        Vector position = {(float) x,(float) y};
        Cell* cell = device_access_grid(grid, position, size);
        switch (cell->type) 
	{
	    case B:
	    {
		lympho_B_action(cell, grid, next_grid, size);
		break;
	    }
	    case Ab:
	    {
		search_antigens(cell, grid, next_grid, size);
	    	break;
	    }
	}
    }
}

__global__
void device_clear_grid(Cell* grid, int subgrid_size, int grid_size, Vector start_position)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < subgrid_size && y < subgrid_size)
    {
        Vector position = {start_position.x + (float) x, start_position.y + (float) y};
        Cell* cell = device_access_grid(grid, position, grid_size);
	    cell->type = FREE;
        cell->position.x = 0;
        cell->velocity.x = 0;
        cell->position.y = 0;
        cell->velocity.y = 0;
        cell->status = INACTIVE;
        for (int i = 0; i < RECEPTOR_SIZE; i++)
        {
            cell->receptor[i] = 0;
        }
    }
}

__global__
void init_cuda_rand(curandState *rand_state, int size)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    if (x < size && y < size)
    {
	long offset = x * size + y;
	curand_init(clock(), 1, 1, &rand_state[offset]);
    }
}


/**
 * Actual simulation function that for n timesteps updates cells in the grid.
 * @param grid grid of the current iteration
 * @param next_grid grid of the next iteration
 * @param options simulation options
 */
__host__
void simulation(Cell* grid, Cell* next_grid, Options options)
{
    dim3 block_size = {16, 16, 1};
    unsigned int grid_length_x = (options.grid_size) / block_size.x + 1;
    unsigned int grid_length_y = (options.grid_size) / block_size.y + 1;
    dim3 number_blocks = {grid_length_x, grid_length_y, 1};
    dim3 thread_per_block = {block_size.x, block_size.y, 1};

    curandState *rand_state;
    cudaMalloc(&rand_state, sizeof(curandState) * options.grid_size * options.grid_size);
    init_cuda_rand<<<number_blocks, thread_per_block>>>(rand_state, options.grid_size);

    for (int t = 0; t < TIMESTEPS; t++)
    {
        timestep<<<number_blocks, thread_per_block>>>(grid, next_grid, rand_state, options.grid_size);
	cudaDeviceSynchronize();
	movement<<<number_blocks, thread_per_block>>>(grid, next_grid, rand_state, options.grid_size);
	cudaDeviceSynchronize();
        /**
         * Swapping grids, so that for every iteration it is used the same grid variable.
         */
	Cell* temp = next_grid;
	next_grid = grid;
	grid = temp;
	Vector start = {0, 0};
	device_clear_grid<<<number_blocks, thread_per_block>>>(next_grid, options.grid_size, options.grid_size, start);
	cudaDeviceSynchronize();
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

    struct timeval t_start_1, t_end_1, t_start_2, t_end_2;

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
    clear_grid(&h_next_grid);
    cudaMemcpy((void*) d_grid.matrix, (void*) h_grid.matrix, malloc_size, cudaMemcpyHostToDevice);
    cudaMemcpy((void*) d_next_grid.matrix, (void*) h_next_grid.matrix, malloc_size, cudaMemcpyHostToDevice);

    //Save cells at the start in a file
    const char* start_file = "start.csv"; 
    save_grid(&h_grid, start_file);

    gettimeofday(&t_start_1, NULL);
    simulation(d_grid.matrix, d_next_grid.matrix, options);
    gettimeofday(&t_end_1, NULL);

    //Save cells in the middle of the simulation, before inserting new antigens.

    cudaMemcpy((void*) h_grid.matrix, (void*) d_grid.matrix, malloc_size, cudaMemcpyDeviceToHost);
    const char* mid_file = "mid.csv";
    save_grid(&h_grid, mid_file);
    insert_antigens(&h_grid);
    cudaMemcpy((void*) d_grid.matrix, (void*) h_grid.matrix, malloc_size, cudaMemcpyHostToDevice);

    gettimeofday(&t_start_2, NULL);
    simulation(d_grid.matrix, d_next_grid.matrix, options);
    gettimeofday(&t_end_2, NULL);

    //Save cells at the end in a file
    cudaMemcpy((void*) h_grid.matrix, (void*) d_grid.matrix, malloc_size, cudaMemcpyDeviceToHost);
    const char* end_file = "end.csv";
    save_grid(&h_grid, end_file);

    printf("T: %lf\n", (WALLTIME(t_end_1) - WALLTIME(t_start_1)) + (WALLTIME(t_end_2) - WALLTIME(t_start_2)));

    free(h_grid.matrix);
    free(h_next_grid.matrix);
    cudaFree(d_grid.matrix);
    cudaFree(d_next_grid.matrix);
    return 0;
}
