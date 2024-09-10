#include <stdio.h>

#include "functions.h"

inline Cell* access_grid(Grid* grid, Vector position)
{
    /**
     * Because of the choice of using a single array grid, accessing it is way much more complicated,
     * but next line has the same meaning of &matrix[position.x][position.y] for a multi-array grid.
     */
    return &grid->matrix[(int)position.x * grid->size + (int)(position).y];
}

void lympho_B_action(Cell* b, Grid* old_grid, Grid* new_grid)
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

void default_action(Cell *cell, Grid *old_grid, Grid *new_grid)
{
    //Empty action for entities with any particular action to do.
}

void search_antigens(Cell* cell, Grid* old_grid, Grid* new_grid)
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
            Cell* other = access_grid(old_grid, current_position);
            if (is_matching_antigen(*cell, *other))
            {
                find_antigen(cell, other);
                return;
            }
        }
    }
}

bool is_matching_antigen(Cell cell, Cell other) 
{
    //Hamming distance of the two receptors has to be greater or equal than the threshold.
    return other.type == Ag && hamming_distance(cell.receptor, other.receptor) >= AFFINITY_MIN;
}

int hamming_distance(char receptor_cell[RECEPTOR_SIZE], char receptor_other[RECEPTOR_SIZE])
{
    int distance = 0;
    char xor[RECEPTOR_SIZE];
    /**
     * Computing xor between every char of the two receptors' arrays.
     */
    for (int i = 0; i < RECEPTOR_SIZE; i++)
    {
        xor[i] = (receptor_cell[i] ^ receptor_other[i]);
    }
    /**
     * Computing the hamming distance as the number of the bit set to 1.
     */
    for (int i = 0; i < RECEPTOR_SIZE; i++) 
    {
        for (int j = 0; j < 8; j++)
        {
            distance += (xor[i] >> j) & 0x1;
        }
    }
    return distance;
}

void correct_position(Vector* position, int size) 
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

void search_lympho_T(Cell* b, Grid* old_grid)
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
            Cell* other = access_grid(old_grid, current_position);
            if (other->type == T)
            {
                b->status = OPERATIVE;
                return;
            }
        }
    }
}

void duplicate(Cell* cell, Grid* old_grid, Grid* new_grid)
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
            Cell* free_cell = access_grid(new_grid, new_position);
            if (free_cell->type == FREE)
            {
                create_duplicate(*cell, free_cell, new_position);
                cell->status = INACTIVE;
                duplicated = true;
            }
        }
    }
}

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

void copy_receptor(unsigned char new_receptor[RECEPTOR_SIZE], unsigned char old_receptor[RECEPTOR_SIZE])
{
    for (int i = 0; i < RECEPTOR_SIZE; i++)
    {
        new_receptor[i] = old_receptor[i];
    }
}

void create_antibodies(Cell* cell, Grid* old_grid, Grid* new_grid)
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
            Cell* free_cell = access_grid(new_grid, new_position);
            if (free_cell->type == FREE)
            {
                create_antibody(*cell, free_cell, new_position);
                created++;
            }
        }
    }
}

void create_antibody(Cell B_cell, Cell* antibody, Vector position)
{
    antibody->action = search_antigens;
    antibody->type = Ab;
    antibody->velocity.x = 0;
    antibody->velocity.y = 0;
    antibody->position = position;
    copy_receptor(antibody->receptor, B_cell.receptor);
}
