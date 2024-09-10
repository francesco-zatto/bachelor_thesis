#include <math.h>
#include <stdlib.h>

#include "physics.h"
#include "functions.h"

/**
 * Function to find a free cell nearby if given starting cell is taken.
 * @param start start cell to search for another free cell
 * @param grid grid where the search is taken
 * @param return free cell near start
 */
static Cell* find_free_cell_nearby(Cell* start, Grid* grid)
{
    for (int s = 1; s <= PROXIMITY_DISTANCE; s++)
    {
        for (int i = -s; i < s; i += s)
        {
            for (int j = -s; j < s; j += s)
            {
                Vector position = {start->position.x + i, start->position.y + j};
                Cell* cell = access_grid(grid, position);
                if (cell->type == FREE)
                {
                    return cell;
                }
            }
        }
    }
    return start;
}

void movement(Cell *cell, Grid *new_grid)
{
    if (cell->type == FREE)
        return;
    
    float box_muller_number[2];
    box_muller(box_muller_number);
    double mass = get_mass(cell->type);
    Vector delta_velocity = 
    {
        .x = langevin_equation(cell->velocity.x, box_muller_number[0], mass),
        .y = langevin_equation(cell->velocity.y, box_muller_number[1], mass)
    };
    cell->velocity.x += delta_velocity.x;
    cell->velocity.y += delta_velocity.y;
    cell->position.x += round(cell->velocity.x * TIMESTEP);
    cell->position.y += round(cell->velocity.y * TIMESTEP);
    correct_position(&(cell->position), new_grid->size);
    Cell* new = access_grid(new_grid, cell->position);
    if (new->type != FREE)
    {
        new = find_free_cell_nearby(new, new_grid);
    }
    *new = *cell;
}

double inline langevin_equation(double velocity, double collision_forces, double mass)
{
    return (-LAMBDA * velocity + collision_forces) / mass * TIMESTEP;
}

double inline get_mass(Type type)
{
    switch (type)
    {
    case B:
    case T:
        return 0.2;
    case Ag:
    case Ab:
        return 0.04;
    }
}

void box_muller(float box_muller_number[2])
{
    float 
        u1 = (float)rand()/(float)(RAND_MAX),
        u2 = (float)rand()/(float)(RAND_MAX);
    box_muller_number[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
    box_muller_number[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2);
}