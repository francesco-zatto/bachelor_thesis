#ifndef _CELL_HEADER
#define _CELL_HEADER

/**
 * Size of a cell's receptor, as a array of chars.
 */
#define RECEPTOR_SIZE 2

/**
 * Minimum requested of affinity in two receptors to make the cells bind, measured in bits.
 */
#define AFFINITY_MIN (RECEPTOR_SIZE * 5)  

/**
 * Maximum distance to limit the region where the cells searches for a free space or antigens.
 */
#define PROXIMITY_DISTANCE 5

/**
 * Number of antigens created by a B lymphocyte when it is operative.
 */
#define NUMBER_CREATED_ANTIGENS 4

#define UNDEFINED -1

/**
 * Grid size used in tests in langevin folder and in tests folder.
 */
#define SIZE 2000

/**
 * Types of cells in the simulation's grid.
 */
typedef enum {B, T, Ag, Ab, FREE} Type;

/**
 * Types of status taken by grid's cells. For now, only B lymphocytes.
 */
typedef enum {INACTIVE, ACTIVE, OPERATIVE} Status;

/**
 * A 2D Vector used for position and velocity.
 */
typedef struct {
    float x;
    float y;
} Vector;

struct Grid;

/**
 * Cell type definition. Every cell has a type, a position, a velocity, a current status, a receptor and an action.
 * Action is a pointer to a function that takes the reference to the cell itself, a pointer to the old grid and a
 * pointer to the next grid.
 */
typedef struct Cell {
    Type type;
    Vector position;
    Vector velocity;
    Status status;
    unsigned char receptor[RECEPTOR_SIZE];
    void (*action)(struct Cell*, struct Grid*, struct Grid*);
} Cell;

/**
 * Grid type definition. A grid contains the matrix of cells and the size of the matrix.
 */
typedef struct Grid {
    Cell* matrix;
    int size;
} Grid;

/**
 * Options type definition. It contains simulation info, like the number of starting cells and the grid size.
 */
typedef struct {
    int total_number_cells;
    int cells_B_number;
    int cells_T_number;
    int ag_number;
    int grid_size;
} Options;

#endif