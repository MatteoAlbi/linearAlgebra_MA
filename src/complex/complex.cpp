#include "complex/complex.hpp"

namespace MA
{

#pragma region constructors

Complex::Complex(double real, double im):
_real(real),
_imag(im)
{}

Complex::Complex(std::pair<double,double> pair):
_real(pair.first),
_imag(pair.second)
{}

#pragma endregion constructors


#pragma region assignment

Complex& Complex::operator=(const Complex& other){
    this->_real = other._real;
    this->_imag = other._imag;
    return *this;
}

Complex& Complex::operator=(Complex const * const other){
    this->_real = other->_real;
    this->_imag = other->_imag;
    return *this;
}

Complex& Complex::operator=(const double& d){
    this->_real = d;
    this->_imag = 0;
    return *this;
}

#pragma endregion assignment


#pragma region set_get

double& Complex::real(){return _real;}
const double& Complex::real() const {return _real;}
void Complex::real(const double & r) {_real = r;}

double& Complex::imag(){return _imag;}
const double& Complex::imag() const {return _imag;}
void Complex::imag(const double & i) {_imag = i;}

#pragma endregion set_get


#pragma region operations

double Complex::magnitude() const {return std::sqrt(_real * _real + _imag * _imag);}
Complex Complex::conj() const {return Complex(_real, -_imag);}

#pragma endregion operations


#pragma region comparators

bool Complex::operator==(const Complex & other) const {return _real == other._real && _imag == other._imag;}
bool Complex::operator==(const Complex * const other) const {return _real == other->_real && _imag == other->_imag;}

bool Complex::operator!=(const Complex & other) const {return _real != other._real || _imag != other._imag;}
bool Complex::operator!=(const Complex * const other) const {return _real != other->_real || _imag != other->_imag;}

#pragma endregion comparators


#pragma region math_operators

void Complex::operator+=(const double & k){this->_real+=k;}
void Complex::operator+=(const Complex & other){
    this->_real += other._real; 
    this->_imag += other._imag; 
}
Complex operator+(const Complex& comp, const double& k){
    return Complex(comp.real() + k, comp.imag());
}
Complex operator+(const double& k, const Complex& comp){
    return Complex(comp.real() + k, comp.imag());
}
Complex operator+(const Complex& comp1, const Complex& comp2){
    return Complex(comp1.real() + comp2.real(), comp1.imag() + comp2.imag());
}


void Complex::operator-=(const double & k) {this->_real-=k;}
void Complex::operator-=(const Complex & other){
    this->_real -= other._real; 
    this->_imag -= other._imag; 
}
Complex operator-(const Complex& comp){
    return Complex(-comp.real(), -comp.imag());
}
Complex operator-(const Complex& comp, const double& k){
    return Complex(comp.real() - k, comp.imag());
}
Complex operator-(const double& k, const Complex& comp){
    return Complex(k - comp.real(), -comp.imag());
}
Complex operator-(const Complex& comp1, const Complex& comp2){
    return Complex(comp1.real() - comp2.real(), comp1.imag() - comp2.imag());
}

void Complex::operator*=(const double & k){
    this->_real *= k; 
    this->_imag *= k; 
}
void Complex::operator*=(const Complex & other){
    double tmp = _real * other._real - _imag * other._imag;
    _imag = _real * other._imag + _imag * other._real;
    _real = tmp;
}
Complex operator*(const Complex& comp, const double& k){
    return Complex(comp.real() * k, comp.imag() * k);
}
Complex operator*(const double& k, const Complex& comp){
    return Complex(comp.real() * k, comp.imag() * k);
}
Complex operator*(const Complex& comp1, const Complex& comp2){
    return Complex(
        comp1.real() * comp2.real() - comp1.imag() * comp2.imag(),
        comp1.real() * comp2.imag() + comp1.imag() * comp2.real()
    );
}

void Complex::operator/=(const double & k){
    if(k == 0.0) throw std::runtime_error("Division by zero");
    this->_real /= k; 
    this->_imag /= k; 
}
void Complex::operator/=(const Complex & other){
    if(other == 0.0) throw std::runtime_error("Division by zero");
    double denominator = other._real * other._real + other._imag * other._imag;
    double tmp = (_real * other._real + _imag * other._imag) / denominator;
    _imag = (_imag * other._real - _real * other._imag) / denominator;
    _real = tmp;
}
Complex operator/(const Complex& comp, const double& k){
    if(k == 0.0) throw std::runtime_error("Division by zero");
    return Complex(comp.real() / k, comp.imag() / k);
}
Complex operator/(const double& k, const Complex& comp){
    if(comp == 0.0) throw std::runtime_error("Division by zero");
    double denominator = comp.real() * comp.real() + comp.imag() * comp.imag();
    return Complex(
        k * comp.real() / denominator,
        - k * comp.imag() / denominator
    );
}
Complex operator/(const Complex& comp1, const Complex& comp2){
    if(comp2 == 0.0) throw std::runtime_error("Division by zero");
    double denominator = comp2.real() * comp2.real() + comp2.imag() * comp2.imag();
    return Complex(
        (comp1.real() * comp2.real() + comp1.imag() * comp2.imag()) / denominator,
        (comp1.imag() * comp2.real() - comp1.real() * comp2.imag()) / denominator
    );
}

#pragma endregion math_operators


std::ostream& operator<<(std::ostream& os, const Complex& comp){
    os << comp.real();
    if(comp.imag() > 0.0) os << "+";
    if(comp.imag() != 0.0) os << comp.imag() << "j";

    return os;
}

std::ostream& operator<<(std::ostream& os, const Complex * comp){
    os << *comp;
    return os;
}

} // namespace MA