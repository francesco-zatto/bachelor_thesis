#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../cell.h"

/**
 * Fluid viscosity. It slows down the cell's movement in the fluid.
 */
#define LAMBDA 0.1

/**
 * Mass of the cell. Depending on it, the cell gains more or less velocity.
 */
#define MASS 2

/**
 * Size of each misuration timestep.
 */
#define TIMESTEP 0.2

/**
 * Size of grid size.
 */
#define SIZE 2000

/**
 * Number of iterations of the Langevin equation simulation.
 */
#define ITERATIONS 1000

#define PI 3.14159265358979323846

/**
 * Default number of cells of the simulation
 */
int cells_number = 2;

/**
 * Using Langevin equation and following a brownian motion inside a fluid, it updates cells coordinates
 * inside 2D grid, so moving with discrete steps.
 * @param cells cells to update their current velocity and position
 * @param n number of cells to update
 */
void update_cell_Langevin(Cell* cells, int n);

/**
 * It generates two random numbers following a normal distribution.
 * @param box_muller_number two numbers following a normal distribution generated by Box-Muller method.
 */
void box_muller(float box_muller_numeber[2]);

/**
 * It save cells positions in an external file.
 * @param positions positions taken by cells inside the grid
 * @param n number of cells in the simulation
 */
void save_positions(Vector* positions, int n);

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    if (argc > 1) 
    {
        cells_number = atoi(argv[1]);
    }
    Cell* cells = calloc(cells_number, sizeof(Cell));
    Vector* positions = malloc(sizeof(Vector) * cells_number * ITERATIONS);
    for (int i = 0; i < ITERATIONS; i++) 
    {
        update_cell_Langevin(cells, cells_number);
        for (int j = 0; j < cells_number; j++) 
        {
            positions[i * cells_number + j] = cells[j].position;
            //printf("%f\t%f\n", cells[j].velocity.x, cells[j].velocity.y);
            printf("%f\t%f\n", positions[i * cells_number + j].x, positions[i * cells_number + j].y);
        }
    }

    save_positions(positions, cells_number);
    free(cells);
    free(positions);
    return 0;
}

void update_cell_Langevin(Cell* cells, int n) 
{
    float box_muller_number[2];
    for (int i = 0; i < n; i++) 
    {
        box_muller(box_muller_number);
        Vector delta_velocity = 
        {
            .x = (-LAMBDA * cells[i].velocity.x + box_muller_number[0]) / MASS * TIMESTEP,
            .y = (-LAMBDA * cells[i].velocity.y + box_muller_number[1]) / MASS * TIMESTEP
        };
        cells[i].velocity.x += delta_velocity.x;
        cells[i].velocity.y += delta_velocity.y;
        cells[i].position.x += round(cells[i].velocity.x * TIMESTEP);
        cells[i].position.y += round(cells[i].velocity.y * TIMESTEP);
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

void save_positions(Vector* positions, int n) 
{
    FILE* file = fopen("data.csv", "w");
    for (int i = 0; i < n; i++) {
        fprintf(file, "x%d;", i);
        fprintf(file, "y%d;", i);
    }
    fprintf(file, "\n");
    for (int i = 0; i < ITERATIONS; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            fprintf(file, "%f;%f;", positions[i * cells_number + j].x, positions[i * cells_number + j].y);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}
