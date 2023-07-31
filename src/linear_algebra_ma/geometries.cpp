#include "linear_algebra_ma/geometries.hpp"

namespace TM
{
    
#pragma region Point

Point::Point(double x, double y): _x(x), _y(y) {}

const double & Point::x() const {return this->_x;}
const double & Point::y() const {return this->_y;}
double & Point::x() {return this->_x;}
double & Point::y() {return this->_y;}
void Point::x(const double & x) {this->_x = x;}
void Point::y(const double & y) {this->_y = y;}

bool Point::operator==(const Point & p) const{
    if(this->isnan() || p.isnan()) return false;
    if(this->_x == p._x && this->_y == p._y) return true;
    else return false;
}

bool Point::operator==(Point const * const p) const{
    return this->operator==(*p);
}

bool Point::operator!=(const Point & p) const{
    return ! this->operator==(p);
}

bool Point::operator!=(Point const * const p) const{
    return ! this->operator==(*p);
}

double Point::distance(const Point & p) const{
    return sqrt(pow(this->_x - p._x, 2) + pow(this->_y - p._y,2));
}

bool Point::isnan() const{
    return (std::isnan(_x) || std::isnan(_y));
}


#pragma endregion Point


#pragma region 

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
    
}

#pragma endregion Line

} // namespace TM
