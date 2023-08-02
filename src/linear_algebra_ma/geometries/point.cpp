#include "linear_algebra_ma/geometries/point.hpp"

namespace MA
{
namespace geometries
{

Point::Point(double x, double y): _x(x), _y(y) {}

const double & Point::x() const {return _x;}
const double & Point::y() const {return _y;}
double & Point::x() {return _x;}
double & Point::y() {return _y;}
std::pair<const double &, const double &> Point::get() const {return {_x, _y};}
void Point::x(const double & x) {_x = x;}
void Point::y(const double & y) {_y = y;}
void Point::set(const double & x, const double & y){ _x = x, _y = y;}

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
    if(this->isnan() || p.isnan()) throw std::invalid_argument("One of the point is NAN");
    if(this->isinf() || p.isinf()) throw std::invalid_argument("One of the point is INF");
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
