#ifndef MA_LINE
#define MA_LINE

#include <utility>
#include <vector>
#include <stdio.h>
#include <string>
#include <math.h>

#include "linear_algebra_ma/matrices.hpp"
#include "linear_algebra_ma/geometries/point.hpp"


namespace MA
{
namespace geometries
{

class Segment{
public:
    Segment(Point p1 = {NAN,NAN}, Point p2 = {NAN,NAN});

    const Point & p1() const;
    const Point & p2() const;
    Point & p1();
    Point & p2();
    std::pair<const Point &,const Point &> get() const;
    void p1(const Point & p1);
    void p2(const Point & p2);
    void set(const Point & p1, const Point & p2);

    bool operator==(const Segment & l) const;
    bool operator==(Segment const * const l) const;

    bool operator!=(const Segment & l) const;
    bool operator!=(Segment const * const l) const;

    bool isnan() const;
    bool isinf() const;

    double length() const;
    double slope() const;

    double distance(const Point & p) const;

    Point intersection(const Segment & l) const;

private:
    Point _p1;
    Point _p2;
};

} // namespace geometries
} // namespace MA

// cout operators
std::ostream& operator<<(std::ostream& os, const MA::geometries::Segment& s);

#endif // MA_LINE