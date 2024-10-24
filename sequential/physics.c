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

#include <math.h>
#include <stdlib.h>
#include <assert.h>

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
    for (int i = -PROXIMITY_DISTANCE; i <= PROXIMITY_DISTANCE; i++)
    {
        for (int j = -PROXIMITY_DISTANCE; j <= PROXIMITY_DISTANCE; j++)
        {
            /**
             * Looking for a nearby position and checking if its free. In that case, that cell will be taken.
             */
            Vector position = {start->position.x + i, start->position.y + j};
            correct_position(&position, grid->size);
            Cell* cell = access_grid(grid, position);
            if (cell->type == FREE)
            {
                return cell;
            }
        }
    }
    return start;
}

void movement(Cell *cell, Grid *new_grid)
{
    //If cell is free, ignore it.
    if (cell->type == FREE)
        return;
    
    /**
     * Computing box muller numbers for forces felt by the body and getting body mass.
     */
    float box_muller_number[2];
    box_muller(box_muller_number);
    double mass = get_mass(cell->type);

    /**
     * Computing deltaV with Langevin equation and then updating cell's position and velocity.
     */
    Vector delta_velocity = 
    {
        .x = langevin_equation(cell->velocity.x, box_muller_number[0], mass),
        .y = langevin_equation(cell->velocity.y, box_muller_number[1], mass)
    };
    cell->velocity.x += delta_velocity.x;
    cell->velocity.y += delta_velocity.y;
    cell->position.x += round(cell->velocity.x * TIMESTEP);
    cell->position.y += round(cell->velocity.y * TIMESTEP);

    /**
     * Checking if the computed position is inside the grid and if it is free.
     */
    correct_position(&(cell->position), new_grid->size);
    Cell* new = access_grid(new_grid, cell->position);
    if (new->type != FREE)
    {
        new = find_free_cell_nearby(new, new_grid);
    }
    new->position = cell->position;
    new->velocity = cell->velocity;
    new->status = cell->status;
    new->type = cell->type;
    new->action = cell->action;
    copy_receptor(new->receptor, cell->receptor);
}

double inline langevin_equation(double velocity, double collision_forces, double mass)
{
    return (-LAMBDA * velocity + collision_forces) / mass * TIMESTEP;
}

double inline get_mass(Type type)
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

void box_muller(float box_muller_number[2])
{
    float 
        u1 = (float)rand()/(float)(RAND_MAX),
        u2 = (float)rand()/(float)(RAND_MAX);
    box_muller_number[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
    box_muller_number[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2);
}