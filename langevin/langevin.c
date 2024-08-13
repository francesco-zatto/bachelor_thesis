#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../cell.h"

/**
 * Fluid viscosity. It slows down the cell's movement in the fluid.
 */
#define LAMBDA 10

/**
 * Mass of the cell. Depending on it, the cell gains more or less velocity.
 */
#define MASS 1

/**
 * Size of each misuration timestep
 */
#define TIMESTEP 0.1

/**
 * Number of iterations of the Langevin equation simulation.
 */
#define ITERATIONS 1000

#define PI 3.14159265358979323846

/**
 * Using Langevin equation and following a brownian motion inside a fluid, it updates cell coordinates.
 * @param cell cell to update its current velocity and position
 */
void update_cell_Langevin(Cell* cell);

/**
 * It generates two random numbers following a normal distribution.
 * @param box_muller_number two numbers following a normal distribution generated by Box-Muller method.
 */
void box_muller(float box_muller_numeber[2]);

/**
 * It save cell positions in an external file.
 */
void save_positions(Vector positions[ITERATIONS]);

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    Vector positions[ITERATIONS];
    Cell cell = {
        .velocity = {
            .x = 0,
            .y = 0
        },
        .position = {
            .x = 0,
            .y = 0
        }
    };
    for (int i = 0; i < ITERATIONS; i++) {
        update_cell_Langevin(&cell);
        positions[i] = cell.position;
    }

    save_positions(positions);
    return 0;
}

void update_cell_Langevin(Cell* cell) 
{
    float box_muller_number[2];
    box_muller(box_muller_number);
    Vector delta_velocity = {
        .x = (-LAMBDA * cell->velocity.x + box_muller_number[0]) / MASS * TIMESTEP,
        .y = (-LAMBDA * cell->velocity.y + box_muller_number[1]) / MASS * TIMESTEP
    };
    cell->velocity.x += delta_velocity.x;
    cell->velocity.y += delta_velocity.y;
    cell->position.x += cell->velocity.x * TIMESTEP;
    cell->position.y += cell->velocity.y * TIMESTEP;
}

void box_muller(float box_muller_number[2])
{
    float 
        u1 = (float)rand()/(float)(RAND_MAX),
        u2 = (float)rand()/(float)(RAND_MAX);
    box_muller_number[0] = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
    box_muller_number[1] = sqrt(-2 * log(u1)) * sin(2 * PI * u2);
}

void save_positions(Vector positions[ITERATIONS]) 
{
    FILE* file = fopen("data.csv", "w");
    for (int i = 0; i < ITERATIONS; i++) {
        fprintf(file, "%f;%f\n", positions[i].x, positions[i].y);
    }
    fclose(file);
}