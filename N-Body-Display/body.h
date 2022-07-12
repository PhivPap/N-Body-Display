#ifndef BODY
#define BODY

#include "util.h"
#include <string>

class Body {
public:
    double mass, vel_x, vel_y;
    Point coords;
    bool ignore = false;

    Body(void) = delete;
    Body(double mass, Point coords, double vel_x, double vel_y);
    void print(void);
};



#endif