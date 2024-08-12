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
    char receptor[RECEPTOR_SIZE];
    void (*action)(void*);
} Cell;