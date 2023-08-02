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

bool Point::isnan() const{
    return (std::isnan(_x) || std::isnan(_y));
}

bool Point::isinf() const{
    return (std::isinf(_x) || std::isinf(_y));
}


double Point::distance(const Point & p) const{
    if(this->isnan() || p.isnan()) throw std::invalid_argument("One of the point is NAN");
    // only one point is INF
    if( (this->isinf() && ! p.isinf()) || 
        (! this->isinf() && p.isinf()) 
    ){
        return INFINITY;
    }
    // both point are infinity
    if(this->isinf() && p.isinf()){
        // one coordinate is infinite and the other same is not
        if( (std::isinf(_x) && ! std::isinf(p._x)) ||
            (std::isinf(_y) && ! std::isinf(p._y)) ||
            (! std::isinf(_x) && std::isinf(p._x)) ||
            (! std::isinf(_y) && std::isinf(p._y))
        ){
            return INFINITY;
        }
        // both coordinate are infinity but opposite sign
        else if( (std::isinf(_x) && std::isinf(p._x) && _x*p._x < 0) ||
                 (std::isinf(_y) && std::isinf(p._y) && _y*p._y < 0)
        ){
            return INFINITY;
        }
        // both coordinate are infinity with same sign
        else if( (std::isinf(_x) && std::isinf(p._x) && _x*p._x > 0) ||
                 (std::isinf(_y) && std::isinf(p._y) && _y*p._y > 0)
        ){
            return NAN;
        }
        // all cases should be covered, but for safety add last else case
        else{
            throw std::runtime_error("Unkown error");
        }
    }

    return sqrt(pow(this->_x - p._x, 2) + pow(this->_y - p._y,2));
}

double Point::slope(const Point & p) const{
    if(this->isnan() || p.isnan()) throw std::invalid_argument("One of the point is NAN");
    if(this->operator==(p)){
        // throw std::runtime_error("The two points are the same: 0/0"); // 0/0 = ?
        return NAN;
    } 
    // both points are infinite
    if(this->isinf() && p.isinf()) {
        // a point has both coordinates = inf
        if( (std::isinf(_x) && std::isinf(_y)) ||
            (std::isinf(p._x) && std::isinf(p._y))
        ){
            // throw std::runtime_error("One point has both coordinates equal to INFINITY: inf/inf"); // inf/inf = ?
            return NAN;
        }
        // points have different coordinates = inf
        else if( (std::isinf(_x) && std::isinf(p._y)) ||
                 (std::isinf(_y) && std::isinf(p._x))
        ){
            // throw std::runtime_error("Points have different inf coordinates: inf/inf"); // inf/inf = ?
            return NAN;
        }
        // both points have x coordinate = inf
        else if( (std::isinf(_x) && std::isinf(p._x))){
            // same sign
            if(_x * p._x > 0){
                // throw std::runtime_error("Points have x coordinates = inf with same sign: inf-inf"); // inf-inf = ?
                return NAN;
            }
            // diff sign
            else{
                return 0; // a/inf = 0
            }
        }
        // both points have y coordinate = inf
        else if( (std::isinf(_y) && std::isinf(p._y))){
            // same sign
            if(_y * p._y > 0){
                // throw std::runtime_error("Points have y coordinates = inf with same sign: inf-inf"); // inf-inf = ?
                return NAN;
            }
            // diff sign
            else{
                return INFINITY; // inf/a = inf
            }
        }
        // all cases should be covered, but for safety add last else case
        else{
            throw std::runtime_error("Unkown error");
        }
    }
    // only this point is infinite
    else if(this->isinf()){
        // both coordinates = inf
        if(std::isinf(_x) && std::isinf(_y)){
            // throw std::runtime_error("This point coordinates are both INFINITY: inf/inf"); // inf/inf = ?
            return NAN;
        }
        // only x coordinate = inf
        else if(std::isinf(_x)){
            return 0; // a/inf = 0
        }
        // only y coordinate = inf
        else if(std::isinf(_y)){
            return INFINITY; // inf/a = inf
        }
        // all cases should be covered, but for safety add last else case
        else{
            throw std::runtime_error("Unkown error");
        }
    } 
    // only second point is infinite
    else if(p.isinf()){
        // both coordinates = inf
        if(std::isinf(p._x) && std::isinf(p._y)){
            // throw std::runtime_error("p point coordinates are both INFINITY: inf/inf"); // inf/inf = ?
            return NAN;
        }
        // only x coordinate = inf
        else if(std::isinf(p._x)){
            return 0; // a/inf = 0
        }
        // only y coordinate = inf
        else if(std::isinf(p._y)){
            return INFINITY; // inf/a = inf
        }
        // all cases should be covered, but for safety add last else case
        else{
            throw std::runtime_error("Unkown error");
        }
    }

    if(_x == p.x()) return INFINITY; // a/0 = inf
    return (_y-p.y()) / (_x-p.x());
}


Matrix Point::to_c_vec() const{
    return Matrix(2,1, {_x,_y});
}

Matrix Point::to_r_vec() const{
    return Matrix(1,2, {_x,_y});
}


Point operator-(const Point & p1, const Point & p2){
    return Point(p1.x() - p2.x(), p1.y() - p2.y());
}

Point operator+(const Point & p1, const Point & p2){
    return Point(p1.x() + p2.x(), p1.y() + p2.y());
}


double distance(const Point & p1, const Point & p2){
    return p1.distance(p2);
}

double slope(const Point & p1, const Point & p2){
    return p1.slope(p2);
}


Matrix c_vec(const Point & p){
    p.to_c_vec();
}

Matrix r_vec(const Point & p){
    return p.to_r_vec();
}


} // namespace geometries
} // namespace MA


std::ostream& operator<<(std::ostream& os, const MA::geometries::Point& p){
    os << "Point(" << p.x() << ", " << p.y() << ")" << std::endl;
}