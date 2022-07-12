#include "body.h"
#include <iostream>


Body::Body(double mass, Point coords, double vel_x, double vel_y):
    mass(mass), coords(coords), vel_x(vel_x), vel_y(vel_y) {}


void Body::print(void){
    std::cout << "mass: " << mass << ", x: " << coords.x << ", y: " << coords.y << ", vel_x: " << vel_x << ", vel_y: " << vel_y << std::endl; 
}