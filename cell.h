#define RECEPTOR_SIZE 2
#define UNDEFINED -1

typedef enum {B, T, Ag, Ab, Free} Type;

typedef enum {Inactive, Active, Operative} Status;

typedef struct {
    float x;
    float y;
} Vector;

typedef struct {
    Type type;
    Vector position;
    Vector velocity;
    Status status;
    unsigned char receptor[RECEPTOR_SIZE];
    void (*action)(void*);
} Cell;

typedef struct {
    int total_number_cells;
    int cells_B_number;
    int cells_T_number;
    int ag_number;
} Options;