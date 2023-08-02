#ifndef MA_POINT
#define MA_POINT

#include <utility>
#include <vector>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdexcept>
#include <exception>

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

    double distance(const Point & p) const;

    bool isnan() const;
    bool isinf() const;

private:
    double _x;
    double _y;

};

} // namespace geometries
} // namespace MA

#endif // MA_POINT