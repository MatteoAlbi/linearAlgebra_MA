#include "linear_algebra_ma/geometries/point.hpp"

namespace MA
{
namespace geometries
{

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
    if(this->isnan() || p.isnan()) return NAN;
    if(this->isinf() || p.isinf()) return INFINITY;
    return sqrt(pow(this->_x - p._x, 2) + pow(this->_y - p._y,2));
}

bool Point::isnan() const{
    return (std::isnan(_x) || std::isnan(_y));
}

bool Point::isinf() const{
    return (std::isinf(_x) || std::isinf(_y));
}

} // namespace geometries
} // namespace MA
