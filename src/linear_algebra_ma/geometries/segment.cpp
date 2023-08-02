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

    // TODO: if projection of point does not belong to the segment
    // take closest extreme
    double seg_length = this->length();
    if(seg_length == 0) return _p2.distance(p);

    return fabs((_p2.x() - _p1.x()) * (_p1.y() - p.y()) - 
                (_p2.y() - _p1.y()) * (_p1.x() - p.x())) 
                / seg_length;
}


bool Segment::point_on_seg(const Point & p) const{
    // A and B are the extremities of the current segment
    // C is the point to check

    // Create the vector AB
    Matrix AB = c_vec(_p2 - _p1);
    // Create the vector AC
    Matrix AC = c_vec(p - _p1);

    // check if points are not aligned using cross prod
    Matrix AB3 = AB | Matrix(1,1, {0});
    Matrix AC3 = AC | Matrix(1,1, {0});
    if(AB3.cross(AC3) != Matrix(3,1,{0,0,0})) return false;

    // Compute the dot product of vectors
    double KAC = AB.dot(AC);
    // Compute the square of the segment length 
    double KAB = AB.dot(AB); 

    if (KAC<0 || KAC>KAB) return false; 
    else return true;
}

Point Segment::intersection(const Segment & l) const{
    // TODO: inf and nan management
    if(this->isnan()) throw std::runtime_error("This segment has NAN values");
    if(this->isinf()) throw std::runtime_error("This segment has INFINITY values");
    if(l.isnan()) throw std::runtime_error("Given segment has NAN values");
    if(l.isinf()) throw std::runtime_error("Given segment has INFINITY values");

    Point q1 = l.p1();
    Point q2 = l.p2();

    Matrix b = c_vec(q1 - _p1);
    Matrix A = c_vec(_p2 - _p1) & c_vec(q1 - q2);

    // non singular A
    if(! A.is_sing()){
        Matrix x = Matrix::solve_ls(A, b);
        if( (0 <= x(0) && x(0) <= 1) && (0 <= x(1) && x(1) <= 1)){
            Matrix tmp = c_vec(_p1) * (1 - x(0)) + c_vec(_p2) * x(0);
            return Point(tmp(0), tmp(1));
        }
        else return Point();
    }
    
    // parallel segments
    //this start inside l
    if(l.point_on_seg(_p1)){
        return _p1;
    }
    //this end inside l
    else if(l.point_on_seg(_p2)){
        if(this->point_on_seg(l._p1)) return l._p1;
        else return l._p2;
    }
    else return Point();
}


std::ostream& operator<<(std::ostream& os, const Segment& s){
    os << "Segment[ (" << s.p1().x() << ", " << s.p1().y() << ") - (" << s.p2().x() << ", " << s.p2().y() << ") ]" << std::endl;
    return os;
}

} // namespace geometries
} // namespace MA