#include "InputOutput.h"

IO::Body::Body(void) :
    id(-1), mass(0), x(0), y(0), vel_x(0), vel_y(0) {}

IO::Body::Body(u32 id, f64 mass, f64 x, f64 y, f64 vel_x, f64 vel_y) :
    id(id), mass(mass), x(x), y(y), vel_x(vel_x), vel_y(vel_y) {}

void IO::Body::copy_from(const IO::Body& b) {
    id = b.id;
    mass = b.mass;
    x = b.x;
    y = b.y;
    vel_x = b.vel_x;
    vel_y = b.vel_y;
}