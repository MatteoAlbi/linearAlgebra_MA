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

    /**
     * @brief check if given point intersect this segment
     * i.e. the point belongs to the segment
     * https://lucidar.me/en/mathematics/check-if-a-point-belongs-on-a-line-segment/
     * @param p given point
     * @return true if point intersect the segment
     */
    bool point_on_seg(const Point & p) const;

    /**
     * @brief check if this segment intersects segment l
     * code from:
     * https://blogs.sas.com/content/iml/2018/07/09/intersection-line-segments.html
     * @param l second segment
     * @return point of intersection
     */
    Point intersection(const Segment & l) const;

private:
    Point _p1;
    Point _p2;
};


// cout operators
std::ostream& operator<<(std::ostream& os, const Segment& s);

} // namespace geometries
} // namespace MA


#endif // MA_LINE