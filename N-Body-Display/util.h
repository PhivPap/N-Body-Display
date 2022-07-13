#ifndef UTIL
#define UTIL
#include "MyTypes.h"

class Point{
public:
    f64 x, y;
    Point(void);
    Point(f64 x, f64 y);
    Point(const Point& p);
    f64 distance_to_2(const Point& p2) const;
    static Point get_center_of_mass(const Point& p1, f64 mass1, const Point& p2, f64 mass2);
};

class Area{
public:
    f64 x1, x2, y1, y2;
    Area(void);
    Area(f64 x1, f64 x2, f64 y1, f64 y2);
    bool contains(const Point& p) const;
    f64 side_length() const;
    f64 diagonal_length_2() const;
    Point get_center() const;
};


#endif