#include "linear_algebra_ma/geometries.hpp"

namespace MA
{

Line::Line(Point p1, Point p2): _p1(p1), _p2(p2) {}

const Point & Line::p1() const {return this->_p1;}
const Point & Line::p2() const {return this->_p2;}
Point & Line::p1() {return this->_p1;}
Point & Line::p2() {return this->_p2;}
void Line::p1(const Point & p1) {this->_p1 = p1;}
void Line::p2(const Point & p2) {this->_p2 = p2;}

bool Line::operator==(const Line & l) const{
    if(this->_p1 == l._p1 && this->_p2 == l._p2) return true;
    else return false;
}
bool Line::operator==(Line const * const l) const{
    return this->operator==(*l);
}

bool Line::operator!=(const Line & l) const{
    return ! this->operator==(l);
}
bool Line::operator!=(Line const * const l) const{
    return ! this->operator==(*l);
}

double Line::slope() const{
    if(_p1.isnan() || _p2.isnan()) return NAN;
    if(_p1 == _p2) return NAN;
    if(_p2.x()-_p1.x() == 0) return INFINITY;
    return (_p2.y()-_p1.y()) / (_p2.x()-_p1.x());
}

double Line::distance(const Point & p) const{
    return fabs((_p2.x() - _p1.x()) * (_p1.y() - p.y()) - 
                (_p2.y() - _p1.y()) * (_p1.x() - p.x())) 
                / sqrt(_p2.distance(_p1));
}

Point Line::intersection(const Line & l) const{
    return l.p1();
}

} // namespace MA
