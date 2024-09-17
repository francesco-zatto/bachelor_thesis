#ifndef _PHYSICS_HEADER
#define _PHYSICS_HEADER

#include <curand.h>
#include <curand_kernel.h>

#include "../cell.h"

/**
 * Fluid viscosity. It slows down the cell's movement in the fluid.
 */
#define LAMBDA 0.1

/**
 * Size of each misuration timestep.
 */
#define TIMESTEP 0.2

#define PI 3.14159

/**
 * Function to compute acceleration using Langevin equation depending on given parameters.
 * @param velocity body current velocity
 * @param collision_forces collision forces felt by the body
 * @param mass body mass
 */
__device__ float langevin_equation(float velocity, float collision_forces, float mass);

/**
 * Getter for body mass depending on its type. For example, lymphocytes are way heavier than antibodies.
 * @param type cell type
 * @return cell's mass
 */
__device__ float get_mass(Type type);

/**
 * Using Langevin equation and following a brownian motion inside a fluid, it updates cells coordinates
 * inside 2D grid, so moving with discrete steps.
 * @param cell cell to update their current velocity and position
 * @param new_grid grid to place cell in next timestep
 * @param rand_state CUDA object that generates random numbers
 */
__device__ void movement(Cell *cell, Grid *new_grid, curandState* rand_state);

/**
 * It generates two random numbers following a normal distribution.
 * @param box_muller_number two numbers following a normal distribution generated by Box-Muller method.
 */
__device__ void box_muller(float box_muller_number[2], curandState* rand_state);

#endif