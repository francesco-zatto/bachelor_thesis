#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cell.h"

#define LAMBDA 0.00001
#define MASS 0.1
#define TIMESTEP 1
#define PI 3.14159265358979323846

float box_muller_number[2];

void update_cell_Langevin(Cell* cell);

void box_muller();

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    Cell cell = {
        .velocity = {
            .x = 0.1,
            .y = 0.1
        },
        .position = {
            .x = 100,
            .y = 100
        }
    };
    for (int i = 0; i < 100; i++) {
        update_cell_Langevin(&cell);
        printf("x: %f \t y: %f \n", cell.position.x, cell.position.y);
    }
    return 0;
}

void update_cell_Langevin(Cell* cell) 
{
    box_muller();
    Vector delta_velocity = {
        .x = (-LAMBDA * cell->velocity.x + box_muller_number[0]) / MASS * TIMESTEP,
        .y = (-LAMBDA * cell->velocity.y + box_muller_number[1]) / MASS * TIMESTEP
    };
    cell->velocity.x += delta_velocity.x;
    cell->velocity.y += delta_velocity.y;
    cell->position.x += cell->velocity.x * TIMESTEP;
    cell->position.y += cell->velocity.y * TIMESTEP;
}

void box_muller()
{
    float box_muller_random[2];
    for (int i = 0; i < 2; i++) {
        box_muller_random[i] = (float)rand()/(float)(RAND_MAX);
    }
    float 
        u1 = box_muller_random[0],
        u2 = box_muller_random[1];
    box_muller_number[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
    box_muller_number[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2);
}
