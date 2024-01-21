#ifndef MA_COMPLEX
#define MA_COMPLEX

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>


namespace MA
{

class Complex{
public:
    Complex(double real = 0, double im = 0);
    Complex(std::pair<double,double> pair);

    Complex& operator=(const Complex& other);
    Complex& operator=(Complex const * const other);
    Complex& operator=(const double& d);

    double& real();
    const double& real() const;
    void real(const double & r);

    double& imag();
    const double& imag() const;
    void imag(const double & i);

    double magnitude() const;
    Complex conj() const;

#pragma region comparators
    bool operator==(const Complex & other) const;
    bool operator==(const Complex * const other) const;

    bool operator!=(const Complex & other) const;
    bool operator!=(const Complex * const other) const;
#pragma endregion comparators

#pragma region math_operators
    void operator+=(const double & k);
    void operator+=(const Complex & other);

    void operator-=(const double & k);
    void operator-=(const Complex & other);

    void operator*=(const double & k);
    void operator*=(const Complex & other);
    
    void operator/=(const double & k);
    void operator/=(const Complex & other);
#pragma endregion math_operators

private:
    double _real;
    double _imag;

};


#pragma region operators
std::ostream& operator<<(std::ostream& os, const Complex& comp);
std::ostream& operator<<(std::ostream& os, const Complex * comp);

Complex operator+(const Complex& comp, const double& k);
Complex operator+(const double& k, const Complex& comp);
Complex operator+(const Complex& comp1, const Complex& comp2);

Complex operator-(const Complex& comp);
Complex operator-(const Complex& comp, const double& k);
Complex operator-(const double& k, const Complex& comp);
Complex operator-(const Complex& comp1, const Complex& comp2);

Complex operator*(const Complex& comp, const double& k);
Complex operator*(const double& k, const Complex& comp);
Complex operator*(const Complex& comp1, const Complex& comp2);

Complex operator/(const Complex& comp, const double& k);
Complex operator/(const double& k, const Complex& comp);
Complex operator/(const Complex& comp1, const Complex& comp2);
#pragma endregion operators


} // namespace MA

#endif //MA_COMPLEX