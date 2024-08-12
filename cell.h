#define RECEPTOR_SIZE 2

enum Type {B, T, Ag, Ab, Free};

enum Status {Inactive, Active, Operative};

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