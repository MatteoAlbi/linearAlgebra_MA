#ifndef MA_GEOMETRIES
#define MA_GEOMETRIES

#include <utility>
#include <vector>
#include <stdio.h>
#include <string>
#include <math.h>

#include "linear_algebra_ma/linearAlgebra.hpp"


namespace MA
{

class Point{
public:
    Point(double x = NAN, double y = NAN);

    const double & x() const;
    const double & y() const;
    double & x();
    double & y();
    void x(const double & x);
    void y(const double & y);

    bool operator==(const Point & p) const;
    bool operator==(Point const * const p) const;

    bool operator!=(const Point & p) const;
    bool operator!=(Point const * const p) const;

    double distance(const Point & p) const;

    bool isnan() const;

private:
    double _x;
    double _y;

};


class Line{
public:
    Line(Point p1 = {NAN,NAN}, Point p2 = {NAN,NAN});

    const Point & p1() const;
    const Point & p2() const;
    Point & p1();
    Point & p2();
    void p1(const Point & p1);
    void p2(const Point & p2);

    bool operator==(const Line & l) const;
    bool operator==(Line const * const l) const;

    bool operator!=(const Line & l) const;
    bool operator!=(Line const * const l) const;

    double slope() const;

    double distance(const Point & p) const;

    Point intersection(const Line & l) const;

private:
    Point _p1;
    Point _p2;
};


} // namespace MA


#endif // MA_GEOMETRIES