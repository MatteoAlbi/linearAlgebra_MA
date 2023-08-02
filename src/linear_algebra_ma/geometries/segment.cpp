#include "linear_algebra_ma/geometries/segment.hpp"

namespace MA
{
namespace geometries
{

Segment::Segment(Point p1, Point p2): _p1(p1), _p2(p2) {}

const Point & Segment::p1() const {return _p1;}
const Point & Segment::p2() const {return _p2;}
Point & Segment::p1() {return _p1;}
Point & Segment::p2() {return _p2;}
std::pair<const Point &,const Point &> Segment::get() const{return {_p1,_p2};}
void Segment::p1(const Point & p1) {_p1 = p1;}
void Segment::p2(const Point & p2) {_p2 = p2;}
void Segment::set(const Point & p1, const Point & p2) {_p1 = p1, _p2 = p2;}

bool Segment::operator==(const Segment & l) const{
    if(_p1 == l._p1 && _p2 == l._p2) return true;
    else return false;
}
bool Segment::operator==(Segment const * const l) const{
    return this->operator==(*l);
}

bool Segment::operator!=(const Segment & l) const{
    return ! this->operator==(l);
}
bool Segment::operator!=(Segment const * const l) const{
    return ! this->operator==(*l);
}

bool Segment::isnan() const{
    return (_p1.isnan() || _p2.isnan());
}

bool Segment::isinf() const{
    return (_p1.isinf() || _p2.isinf());
}

double Segment::length() const{
    return _p1.distance(_p2);
}

double Segment::slope() const{
    return _p1.slope(_p2);
}

double Segment::distance(const Point & p) const{
    // TODO: inf and nan management
    if(this->isnan()) throw std::runtime_error("This segment has NAN values");
    if(this->isinf()) throw std::runtime_error("This segment has INFINITY values");

    double line_points_dist = _p2.distance(_p1);
    if(line_points_dist == 0) return _p2.distance(p);

    return fabs((_p2.x() - _p1.x()) * (_p1.y() - p.y()) - 
                (_p2.y() - _p1.y()) * (_p1.x() - p.x())) 
                / sqrt(line_points_dist);
}

Point Segment::intersection(const Segment & l) const{
    // TODO: inf and nan management
    if(this->isnan()) throw std::runtime_error("This segment has NAN values");
    if(this->isinf()) throw std::runtime_error("This segment has INFINITY values");
    
    return l.p1();
}

} // namespace geometries
} // namespace MA
