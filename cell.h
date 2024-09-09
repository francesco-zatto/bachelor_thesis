#ifndef _CELL_HEADER
#define _CELL_HEADER

#define RECEPTOR_SIZE 2
#define AFFINITY_MIN (RECEPTOR_SIZE * 6)  
#define PROXIMITY_DISTANCE 5
#define NUMBER_CREATED_ANTIGENS 8
#define UNDEFINED -1
#define SIZE 2000

typedef enum {B, T, Ag, Ab, FREE} Type;

typedef enum {INACTIVE, ACTIVE, OPERATIVE} Status;

typedef struct {
    float x;
    float y;
} Vector;

struct Grid;

typedef struct Cell {
    Type type;
    Vector position;
    Vector velocity;
    Status status;
    unsigned char receptor[RECEPTOR_SIZE];
    void (*action)(struct Cell*, struct Grid*, struct Grid*);
} Cell;

typedef struct Grid {
    Cell* matrix;
    int size;
} Grid;

typedef struct {
    int total_number_cells;
    int cells_B_number;
    int cells_T_number;
    int ag_number;
    int grid_size;
} Options;

#endif