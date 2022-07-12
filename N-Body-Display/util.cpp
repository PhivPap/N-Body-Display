#include "util.h"
#include <cmath>



/* Point methods */

Point::Point(void) {}

Point::Point(double x, double y) : x(x), y(y) {}

Point::Point(const Point& p){
    x = p.x;
    y = p.y;
}

double Point::distance_to_2(const Point& p) const {
    return (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y);
}

Point Point::get_center_of_mass(const Point& p1, double mass1, const Point& p2, double mass2){
    double total_mass = mass1 + mass2;
    if (total_mass == 0)
        return Point(0, 0);
    return Point((p1.x * mass1 + p2.x * mass2) / total_mass, (p1.y * mass1 + p2.y * mass2) / total_mass);
}

/* Area methods */

Area::Area(void) {}

Area::Area(double x1, double x2, double y1, double y2) :
    x1(x1), x2(x2), y1(y1), y2(y2) {}

bool Area::contains(const Point& p) const {
    return p.x >= x1 && p.x <= x2 && p.y >= y1 && p.y <= y2;
}

double Area::side_length() const {
    return x2 - x1;
}

double Area::diagonal_length_2() const {
    return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
}

Point Area::get_center() const {
    return Point((x1 + x2) / 2, (y1 + y2) / 2);
}

