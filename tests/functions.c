#include "../cell.h"
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char const *argv[])
{
    Cell* b = malloc(sizeof(Cell));
    b->action = lympho_B_action;
    b->type = B;
    b->receptor[0] = 0x0f, b->receptor[1] = 0xf;
    printf("%x %x\n", b->receptor[0] & 0xff, b->receptor[1] & 0xff);
    Cell* ag = malloc(sizeof(Cell));
    ag->type = Ag;
    ag->receptor[0] = 0xff, ag->receptor[1] = 0xff;
    printf("%x %x\n", ag->receptor[0] & 0xff, ag->receptor[1] & 0xff);
    printf("%d\n", hamming_distance(b->receptor, ag->receptor));
    printf("%d\n", is_matching_antigen(*b, *ag));
    free(b);
    free(ag);
    return 0;
}


void lympho_B_action(Cell* b, Cell* grid)
{
    switch (b->status)
    {
        case INACTIVE: 
        {
            search_antigens(b, grid);
            break;
        }
        case ACTIVE: 
        {
            break;
        }
        case OPERATIVE:
        {
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
                .y = cell->position.x + j 
            };
            correct_position(&current_position);
            Cell* other = &grid[(int)current_position.x * SIZE + (int)(current_position).y];
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
    printf("%x %x\n", xor[0] && 0xff, xor[1] & 0xff);
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
        position->x = SIZE - position->x;
    }
    if (position->y >= SIZE)
    {
        position->y -= SIZE;
    }
    if (position->y < 0)
    {
        position->y = SIZE - position->y;
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

