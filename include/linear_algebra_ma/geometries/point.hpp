#ifndef MA_POINT
#define MA_POINT

#include <utility>
#include <vector>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdexcept>
#include <exception>

#include "linear_algebra_ma/matrices.hpp"

namespace MA
{
namespace geometries
{

class Point{
public:
    Point(double x = NAN, double y = NAN);

    const double & x() const;
    const double & y() const;
    double & x();
    double & y();
    std::pair<const double &, const double &> get() const;
    void x(const double & x);
    void y(const double & y);
    void set(const double & x, const double & y);

    bool operator==(const Point & p) const;
    bool operator==(Point const * const p) const;

    bool operator!=(const Point & p) const;
    bool operator!=(Point const * const p) const;

    bool isnan() const;
    bool isinf() const;

    double distance(const Point & p) const;
    double slope(const Point & p) const;

    Matrix to_c_vec() const;
    Matrix to_r_vec() const;

private:
    double _x;
    double _y;

};

Point operator-(const Point & p1, const Point & p2);
Point operator+(const Point & p1, const Point & p2);

double distance(const Point & p1, const Point & p2);
double slope(const Point & p1, const Point & p2);

Matrix c_vec(const Point & p);
Matrix r_vec(const Point & p);


} // namespace geometries
} // namespace MA

// cout operators
std::ostream& operator<<(std::ostream& os, const MA::geometries::Point& p);

#endif // MA_POINT