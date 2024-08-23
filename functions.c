#include "cell.h"
#include "functions.h"

#include <stdio.h>

inline Cell* access_grid(Cell* grid, Vector position)
{
    return &grid[(int)position.x * SIZE + (int)(position).y];
}

void lympho_B_action(Cell* b, Cell* old_grid, Cell* new_grid)
{
    switch (b->status)
    {
        case INACTIVE: 
        {
            search_antigens(b, old_grid);
            break;
        }
        case OPERATIVE:
        {
            duplicate(b, old_grid, new_grid);
            create_antibodies(b, old_grid, new_grid);
            break;
        }
    }
}

void search_antigens(Cell* cell, Cell* grid)
{
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE; j++) 
        {
            Vector current_position = 
            {
                .x = cell->position.x + i,
                .y = cell->position.y + j 
            };
            correct_position(&current_position);
            Cell* other = access_grid(grid, current_position);
            if (is_matching_antigen(*cell, *other))
            {
                find_antigen(cell, other);
            }
        }
    }
}

bool is_matching_antigen(Cell cell, Cell other) 
{
    return other.type == Ag && hamming_distance(cell.receptor, other.receptor) >= AFFINITY_MIN;
}

int hamming_distance(char receptor_cell[RECEPTOR_SIZE], char receptor_other[RECEPTOR_SIZE])
{
    int distance = 0;
    char xor[RECEPTOR_SIZE];
    for (int i = 0; i < RECEPTOR_SIZE; i++)
    {
        xor[i] = (receptor_cell[i] ^ receptor_other[i]);
    }
    for (int i = 0; i < RECEPTOR_SIZE; i++) 
    {
        for (int j = 0; j < 8; j++)
        {
            distance += (xor[i] >> j) & 0x1;
        }
    }
    return distance;
}

void correct_position(Vector* position) {
    if (position->x >= SIZE)
    {
        position->x -= SIZE;
    }
    if (position->x < 0)
    {
        position->x = SIZE + position->x;
    }
    if (position->y >= SIZE)
    {
        position->y -= SIZE;
    }
    if (position->y < 0)
    {
        position->y = SIZE + position->y;
    }
}

void find_antigen(Cell* cell, Cell* other)
{
    switch (cell->type)
    {
        case B:
        {
            cell->status = ACTIVE;
            break;
        }
        case Ag:
        {
            other->type = FREE;
            break;
        }
    }
}

void duplicate(Cell* cell, Cell* old_grid, Cell* new_grid)
{
    bool duplicated = false;
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE && !duplicated; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE && !duplicated; j++) 
        {

            Vector new_position = 
            {
                .x = cell->position.x + i,
                .y = cell->position.y + j
            };
            correct_position(&new_position);
            Cell* free_cell = access_grid(old_grid, new_position);
            if (free_cell->type == FREE)
            {
                create_duplicate(*cell, free_cell, new_position);
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

void create_antibodies(Cell* cell, Cell* old_grid, Cell* new_grid)
{

}


