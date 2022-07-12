#include "InputOutput.h"

IO::Body::Body(void) :
    id(-1), mass(0), x(0), y(0), vel_x(0), vel_y(0) {}

IO::Body::Body(uint32_t id, double mass, double x, double y, double vel_x, double vel_y) :
    id(id), mass(mass), x(x), y(y), vel_x(vel_x), vel_y(vel_y) {}

void IO::Body::copy_from(const IO::Body& b) {
    id = b.id;
    mass = b.mass;
    x = b.x;
    y = b.y;
    vel_x = b.vel_x;
    vel_y = b.vel_y;
}