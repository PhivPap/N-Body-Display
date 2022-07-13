#ifndef BODY
#define BODY

#include "util.h"
#include <string>
#include "MyTypes.h"

class Body {
public:
    f64 mass, vel_x, vel_y;
    Point coords;
    bool ignore = false;

    Body(void) = delete;
    Body(f64 mass, Point coords, f64 vel_x, f64 vel_y);
    void print(void);
};



#endif