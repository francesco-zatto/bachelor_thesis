#include <math.h>
#include <stdlib.h>

#include "physics.h"
#include "functions.h"

void movement(Cell *cell, Cell *new_grid)
{
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
    correct_position(&(cell->position));
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
        return 1;
    case Ag:
    case Ab:
        return 0.08;
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