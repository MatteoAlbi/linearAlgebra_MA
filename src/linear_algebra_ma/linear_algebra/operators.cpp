#include "linear_algebra_ma/linearAlgebra.hpp"

namespace MA
{

#pragma region access
double& Matrix::operator()(const uint & r, const uint & c){
    if(r >= this->_r) throw out_of_range("Row index out of range");
    if(c >= this->_c) throw out_of_range("Col index out of range");
    return this->_v[r * this->_c + c];
}

const double& Matrix::operator()(const uint & r, const uint & c) const{
    if(r >= this->_r) throw out_of_range("Row index out of range");
    if(c >= this->_c) throw out_of_range("Col index out of range");
    return this->_v[r * this->_c + c];
}

double& Matrix::operator()(const uint & i){
    if(this->_c == 1){
        if(i >= this->_r) throw out_of_range("Index out of range");
        return this->operator()(i,0);
    } 
    else if(this->_r == 1){
        if(i >= this->_c) throw out_of_range("Index out of range");
        return this->operator()(0,i);
    }
    else{
        if(i >= this->_r || i >= this->_c) throw out_of_range("Index out of range");
        return this->operator()(i,i);
    }
}

const double& Matrix::operator()(const uint & i) const{
    if(this->_c == 1){
        if(i >= this->_r) throw out_of_range("Index out of range");
        return this->operator()(i,0);
    } 
    else if(this->_r == 1){
        if(i >= this->_c) throw out_of_range("Index out of range");
        return this->operator()(0,i);
    }
    else{
        if(i >= this->_r || i >= this->_c) throw out_of_range("Index out of range");
        return this->operator()(i,i);
    }
}

Matrix Matrix::operator()(pair<uint,uint> rs, const uint & c) const{
    if(rs.second >= this->_r) throw out_of_range("Row index out of range");
    if(c >= this->_c) throw out_of_range("Col index out of range");
    if(rs.first > rs.second) throw invalid_argument("Row first element must be <= of second");

    Matrix ret = Matrix(rs.second - rs.first + 1, 1);

    for(uint i=0; i<ret._r; i++){
        ret(i,0) = this->operator()(i + rs.first, c);
    }

    return ret;
}

Matrix Matrix::operator()(const uint & r, pair<uint,uint> cs) const{
    if(r >= this->_r) throw out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw out_of_range("Col index out of range");
    if(cs.first > cs.second) throw invalid_argument("Col first element must be <= of second");

    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;
    Matrix ret = Matrix(1, cs.second - cs.first + 1);

    for(uint j=0; j<ret._c; j++){
        ret(0,j) = this->operator()(r, j + cs.first);
    }

    return ret;
}

Matrix Matrix::operator()(pair<uint,uint> rs, pair<uint,uint> cs) const{
    if(rs.second >= this->_r) throw out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw out_of_range("Col index out of range");
    if(rs.first > rs.second) throw invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw invalid_argument("Col first element must be <= of second");

    Matrix ret = Matrix(rs.second - rs.first + 1, cs.second - cs.first + 1);

    for(uint i=0; i<ret._r; i++){
        for(uint j=0; j<ret._c; j++){
            ret(i,j) = this->operator()(i + rs.first, j + cs.first);
        }
    }

    return ret;
}
#pragma endregion access

#pragma region assign
Matrix& Matrix::operator=(const Matrix& m){
    if (this == &m) return *this;

    // set dimensions
    uint prev_size = this->size();
    this->_r = m._r;
    this->_c = m._c;

    // redimension
    if(prev_size != this->size()){
        if(this->_v != nullptr) delete[] this->_v;
        this->_v = new double[this->size()];
    }
    // set values
    std::copy(m._v, m._v + m.size(), this->_v);

    return *this;
}

Matrix& Matrix::operator=(Matrix const * const m){
    return this->operator=(*m);
}
#pragma endregion assign

#pragma region comparators
bool Matrix::operator==(const Matrix & m) const{
    // check shape
    if(this->getR() != m.getR()) return false;
    if(this->getC() != m.getC()) return false;
    // check values
    for(uint i=0; i<m.getR(); i++) for(uint j=0; j<m.getC(); j++){
        if((*this)(i,j) != m(i,j)) return false;
    }

    return true;
}

bool Matrix::operator==(const Matrix * const m) const{
    return this->operator==(*m);
}

bool Matrix::operator!=(const Matrix & m) const{
    return ! this->operator==(m);
}

bool Matrix::operator!=(const Matrix * const m) const{
    return ! this->operator==(*m);
}
#pragma endregion comparators

#pragma region sum
void Matrix::operator+=(const double & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) += k;
        }
    }
}

void Matrix::operator+=(const int & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) += k;
        }
    }
}

void Matrix::operator+=(const Matrix & m){
    if(this->_r != m._r || this->_c != m._c) throw invalid_argument("Matrices' shapes don't match");
    
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) += m(i,j);
        }
    }
}

Matrix operator+(const Matrix& m, const double& k){
    Matrix ret = m;
    ret+=k;

    return ret;
}

Matrix operator+(const Matrix& m, const int& k){
    Matrix ret = m;
    ret+=k;

    return ret;
}

Matrix operator+(const Matrix& m1, const Matrix& m2){
    if(m1.getR() != m2.getR() || m1.getC() != m2.getC()) throw invalid_argument("Matrices' shapes don't match");
    
    Matrix ret = m1;
    ret+=m2;

    return ret;
}
#pragma endregion sum

#pragma region subtract
void Matrix::operator-=(const double & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) -= k;
        }
    }
}

void Matrix::operator-=(const int & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) -= k;
        }
    }
}

void Matrix::operator-=(const Matrix & m){
    if(this->_r != m._r || this->_c != m._c) throw invalid_argument("Matrices' shapes don't match");
    
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) -= m(i,j);
        }
    }
}

Matrix operator-(const Matrix& m, const double& k){
    Matrix ret = m;
    ret-=k;

    return ret;
}

Matrix operator-(const Matrix& m, const int& k){
    Matrix ret = m;
    ret-=k;

    return ret;
}

Matrix operator-(const Matrix& m1, const Matrix& m2){
    if(m1.getR() != m2.getR() || m1.getC() != m2.getC()) throw invalid_argument("Matrices' shapes don't match");
    
    Matrix ret = m1;
    ret-=m2;

    return ret;
}
#pragma endregion subtract

#pragma region multiply
void Matrix::operator*=(const double & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) *= k;
        }
    }
}

void Matrix::operator*=(const int & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) *= k;
        }
    }
}

void Matrix::operator*=(const Matrix & m){
    this->operator=((*this) * m);
}

Matrix operator*(const Matrix& m, const double& k){
    Matrix ret = m;
    ret*=k;

    return ret;
}

Matrix operator*(const Matrix& m, const int& k){
    Matrix ret = m;
    ret*=k;

    return ret;
}

Matrix operator*(const Matrix& m1, const Matrix& m2){
    if(m1.getC() != m2.getR()) throw invalid_argument("Matrices' shapes don't match");
    
    Matrix ret = Matrix(m1.getR(), m2.getC());

    for (uint i = 0; i<m1.getR(); i++) {
        for (uint j = 0; j<m2.getC(); j++) {
            ret(i,j) = 0;
            for (uint k = 0; k<m1.getC(); k++) {
                ret(i,j) += m1(i,k) * m2(k,j);
            }
        }
    }

    return ret;
}
#pragma endregion multiply

#pragma region right_divide
void Matrix::operator/=(const double & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) /= k;
        }
    }
}

void Matrix::operator/=(const int & k){
    for(uint i=0; i<this->_r; i++){
        for(uint j=0; j<this->_c; j++){
            this->operator()(i,j) /= k;
        }
    }
}

/*
void Matrix::operator/=(const Matrix & m){
    Matrix inv;
    try
    {
        if(m.getC() == m.getR()) inv = m.inv();
        else inv = m.pinv_left();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        throw invalid_argument("Matrix m not invertible");
    }    

    this->operator*=(inv);
}
*/

Matrix operator/(const Matrix& m, const double& k){
    Matrix ret = m;
    ret/=k;

    return ret;
}

Matrix operator/(const Matrix& m, const int& k){
    Matrix ret = m;
    ret/=k;

    return ret;
}

/*
Matrix operator/(const Matrix& m1, const Matrix& m2){
    Matrix inv;
    try
    {
        if(m2.getC() == m2.getR()) inv = m2.inv();
        else inv = m2.pinv_left();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        throw invalid_argument("Matrix m2 not invertible");
    }    

    return m1*inv;
}
*/
#pragma endregion right_divide

#pragma region concatenate
void Matrix::operator&=(const Matrix & m){
    this->operator=((*this) & m);
}

void Matrix::operator|=(const Matrix & m){
    this->operator=((*this) | m);
}

Matrix operator&(const Matrix& m1, const Matrix& m2){
    if(m1.getR() != m2.getR()) throw invalid_argument("Matrices must have same number of rows");

    Matrix ret = Matrix(m1.getR(), m1.getC()+m2.getC());

    for(uint i=0; i<ret.getR(); i++){
        for(uint j=0; j<ret.getC(); j++){
            if(j < m1.getC()) ret(i,j) = m1(i,j);
            else ret(i,j) = m2(i,j-m1.getC());
        }
    }

    return ret;
}

Matrix operator|(const Matrix& m1, const Matrix& m2){
    if(m1.getC() != m2.getC()) throw invalid_argument("Matrices must have same number of columns");

    Matrix ret = Matrix(m1.getR()+m2.getR(), m1.getC());

    for(uint i=0; i<ret.getR(); i++){
        for(uint j=0; j<ret.getC(); j++){
            if(i < m1.getR()) ret(i,j) = m1(i,j);
            else ret(i,j) = m2(i-m1.getR(),j);
        }
    }

    return ret;
}
#pragma endregion concatenate


ostream& operator<<(ostream& os, const Matrix& m){
    os << "Matrix(" << m.getR() << "x" << m.getC() << ")" << endl;
    for(uint i=0; i<m.getR(); i++){
        if(i>0) os << endl;
        for(uint j=0; j<m.getC(); j++){
            if(j>0) os << " ";
            os << m(i,j);
        }
    }
    
    return os;
}

ostream& operator<<(ostream& os, const Matrix * m){
    os << "Matrix(" << m->getR() << "x" << m->getC() << ")" << endl;
    for(uint i=0; i<m->getR(); i++){
        if(i>0) os << endl;
        for(uint j=0; j<m->getC(); j++){
            if(j>0) os << " ";
            os << m->operator()(i,j);
        }
    }
    
    return os;
}



} // namespace MA


