#ifndef MA_MATRICES_OP_TPP
#define MA_MATRICES_OP_TPP

#include "linear_algebra_ma/matrices.hpp"

namespace MA
{

#pragma region access

template<typename T>
T& Matrix<T>::operator()(const uint & r, const uint & c){
    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(c >= this->_c) throw std::out_of_range("Col index out of range");
    return this->_v[r * this->_c + c];
}

template<typename T>
const T& Matrix<T>::operator()(const uint & r, const uint & c) const{
    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(c >= this->_c) throw std::out_of_range("Col index out of range");
    return this->_v[r * this->_c + c];
}

template<typename T>
T& Matrix<T>::operator()(const uint & i){
    if(this->_c == 1){
        if(i >= this->_r) throw std::out_of_range("Index out of range");
        return this->operator()(i,0);
    } 
    else if(this->_r == 1){
        if(i >= this->_c) throw std::out_of_range("Index out of range");
        return this->operator()(0,i);
    }
    else{
        if(i >= this->_r || i >= this->_c) throw std::out_of_range("Index out of range");
        return this->operator()(i,i);
    }
}

template<typename T>
const T& Matrix<T>::operator()(const uint & i) const{
    if(this->_c == 1){
        if(i >= this->_r) throw std::out_of_range("Index out of range");
        return this->operator()(i,0);
    } 
    else if(this->_r == 1){
        if(i >= this->_c) throw std::out_of_range("Index out of range");
        return this->operator()(0,i);
    }
    else{
        if(i >= this->_r || i >= this->_c) throw std::out_of_range("Index out of range");
        return this->operator()(i,i);
    }
}

template<typename T>
Matrix<T> Matrix<T>::operator()(uu_pair rs, const uint & c) const{
    if(rs.second >= this->_r) throw std::out_of_range("Row index out of range");
    if(c >= this->_c) throw std::out_of_range("Col index out of range");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second");

    // allow to use {} for full row selection
    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;

    Matrix ret = Matrix(rs.second - rs.first + 1, 1);

    for(uint i=0; i<ret._r; ++i){
        ret(i,0) = this->operator()(i + rs.first, c);
    }

    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::operator()(const uint & r, uu_pair cs) const{
    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw std::out_of_range("Col index out of range");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    // allow to use {} for full col selection
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    Matrix ret = Matrix(1, cs.second - cs.first + 1);

    for(uint j=0; j<ret._c; ++j){
        ret(0,j) = this->operator()(r, j + cs.first);
    }

    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::operator()(uu_pair rs, uu_pair cs) const{
    if(rs.second >= this->_r) throw std::out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw std::out_of_range("Col index out of range");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    // allow to use {} for full row/col selection
    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    Matrix ret = Matrix(rs.second - rs.first + 1, cs.second - cs.first + 1);

    for(uint i=0; i<ret._r; ++i){
        for(uint j=0; j<ret._c; ++j){
            ret(i,j) = this->operator()(i + rs.first, j + cs.first);
        }
    }

    return ret;
}

#pragma endregion access


#pragma region assign

template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> m){
    swap(*this, m);
    return *this;
}

#pragma endregion assign


#pragma region comparators

template<typename U, typename V>
std::enable_if_t<!is_complex<U>::value && !is_complex<V>::value, bool>
operator==(const Matrix<U> & m1, const Matrix<V> & m2) {
    // check shape
    if(m1.r() != m2.r()) return false;
    if(m1.c() != m2.c()) return false;
    // check values
    for(uint i=0; i<m1.r(); ++i) for(uint j=0; j<m1.c(); ++j){
        double err = abs(m1(i,j) - m2(i,j));
        if( err > Matrix<U>::get_epsilon()) return false;
    }
    return true;
}

template<typename U, typename V>
std::enable_if_t<is_complex<U>::value || is_complex<V>::value, bool>
operator==(const Matrix<U> & m1, const Matrix<V> & m2) {
    // check shape
    if(m1.r() != m2.r()) return false;
    if(m1.c() != m2.c()) return false;
    // check values
    for(uint i=0; i<m1.r(); ++i) for(uint j=0; j<m1.c(); ++j){
        std::complex<double> tmp1 = m1(i,j);
        std::complex<double> tmp2 = m2(i,j);
        double err_real = abs(tmp1.real() - tmp2.real());
        double err_imag = abs(tmp1.imag() - tmp2.imag());
        if( err_real > Matrix<U>::get_epsilon() || err_imag > Matrix<U>::get_epsilon()){ 
            using namespace std;
            cout << "values in " << uu_pair{i,j} << " are different: " << endl
                 << "m1: " << m1(i,j) 
                 << ", m2: " << m2(i,j) 
                 << ", diff: (" << err_real << ", " << err_imag << ")" << endl;
            return false;
        }
    }
    return true;
}

template<typename U, typename V>
bool operator!=(const Matrix<U> & m1, const Matrix<V> & m2){
    return ! (m1==m2);
}

#pragma endregion comparators


#pragma region sum

template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator+=(const V & k){
    for(uint i=0; i<this->size(); ++i){
        this->_v[i] += k;
    }
    return *this;
}

template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator+=(const Matrix<V> & m){
    if(this->_r != m.r() || this->_c != m.c()) throw std::invalid_argument("Matrices' shapes don't match");

    for(uint i=0; i<this->size(); ++i){
        this->_v[i] += m.v()[i];
    }
    return *this;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator+(const Matrix<U>& m, const V& k){
    Matrix<double> ret = m;
    ret+=k;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator+(const Matrix<U>& m, const V& k){
    Matrix<std::complex<double>> ret = m;
    ret+=k;
    return ret;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator+(const U& k, const Matrix<V>& m){
    Matrix<double> ret = m;
    ret+=k;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator+(const U& k, const Matrix<V>& m){
    Matrix<std::complex<double>> ret = m;
    ret+=k;
    return ret;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator+(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.r() != m2.r() || m1.c() != m2.c()) throw std::invalid_argument("Matrices' shapes don't match");
    Matrix<double> ret = m1;
    ret+=m2;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator+(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.r() != m2.r() || m1.c() != m2.c()) throw std::invalid_argument("Matrices' shapes don't match");
    Matrix<std::complex<double>> ret = m1;
    ret+=m2;
    return ret;
}

#pragma endregion sum


#pragma region subtract

template<typename T>
template<typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator-=(const V & k){
    for(uint i=0; i<this->_r; ++i){
        for(uint j=0; j<this->_c; ++j){
            this->operator()(i,j) -= k;
        }
    }
    return *this;
}

template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator-=(const Matrix<V> & m){
    if(this->_r != m.r() || this->_c != m.c()) throw std::invalid_argument("Matrices' shapes don't match");

    for(uint i=0; i<this->size(); ++i){
        this->_v[i] -= m.v()[i];
    }
    return *this;
}


template<typename U>
Matrix<U> operator-(const Matrix<U>& m){
    Matrix ret = m;
    for(uint i=0; i<ret.size(); ++i){
        ret._v[i] = -ret._v[i];
    }
    return ret;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator-(const Matrix<U>& m, const V& k){
    Matrix<double> ret = m;
    ret-=k;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator-(const Matrix<U>& m, const V& k){
    Matrix<std::complex<double>> ret = m;
    ret-=k;
    return ret;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator-(const U& k, const Matrix<V>& m){
    return -m + k;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator-(const U& k, const Matrix<V>& m){
    return -m + k;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator-(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.r() != m2.r() || m1.c() != m2.c()) throw std::invalid_argument("Matrices' shapes don't match");
    Matrix<double> ret = m1;
    ret-=m2;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator-(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.r() != m2.r() || m1.c() != m2.c()) throw std::invalid_argument("Matrices' shapes don't match");
    Matrix<std::complex<double>> ret = m1;
    ret-=m2;
    return ret;
}

#pragma endregion subtract


#pragma region multiply

template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator*=(const V & k){
    for(uint i=0; i<this->size(); ++i){
        this->_v[i] *= k;
    }
    return *this;
}

template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator*=(const Matrix<V> & m){
    this->operator=((*this) * m);
    return *this;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator*(const Matrix<U>& m, const V& k){
    Matrix<double> ret = m;
    ret*=k;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator*(const Matrix<U>& m, const V& k){
    Matrix<std::complex<double>> ret = m;
    ret*=k;
    return ret;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator*(const U& k, const Matrix<V>& m){
    Matrix<double> ret = m;
    ret*=k;
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator*(const U& k, const Matrix<V>& m){
    Matrix<std::complex<double>> ret = m;
    ret*=k;
    return ret;
}


template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator*(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.c() != m2.r()) throw std::invalid_argument("Matrices' shapes don't match");
    
    Matrix<double> ret(m1.r(), m2.c());

    for (uint i = 0; i<m1.r(); ++i) {
        for (uint j = 0; j<m2.c(); ++j) {
            ret(i,j) = 0.0;
            for (uint k = 0; k<m1.c(); ++k) {
                ret(i,j) += m1(i,k) * m2(k,j);
            }
        }
    }

    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator*(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.c() != m2.r()) throw std::invalid_argument("Matrices' shapes don't match");
    
    Matrix<std::complex<double>> ret(m1.r(), m2.c());

    for (uint i = 0; i<m1.r(); ++i) {
        for (uint j = 0; j<m2.c(); ++j) {
            ret(i,j) = 0.0;
            for (uint k = 0; k<m1.c(); ++k) {
                ret(i,j) += m1(i,k) * m2(k,j);
            }
        }
    }

    return ret;
}

#pragma endregion multiply

// todo
#pragma region divide

template<typename T>
Matrix<T>& Matrix<T>::operator/=(const double & k){
    for(uint i=0; i<this->_r; ++i){
        for(uint j=0; j<this->_c; ++j){
            this->operator()(i,j) /= k;
        }
    }
    return *this;
}

// Matrix& Matrix<T>::operator/=(const Matrix & m){
//    this->operator=((*this) / m);
//    return *this;
// }

template<typename T>
Matrix<T> operator/(const Matrix<T>& m, const double& k){
    Matrix<T> ret(m);
    ret/=k;

    return ret;
}

// Matrix operator/(const double& k, const Matrix& m){
//     if(m.r() == m.c()) return m.inv() * k;
//     else if(m.r() > m.c()) return m.pinv_left() * k;
//     else return m.pinv_right() * k;
// }

// Matrix operator/(const Matrix& m1, const Matrix& m2){
//     return Matrix<T>::matrix_r_divide(m1,m2);
// }

#pragma endregion divide


#pragma region concatenate

template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator&=(const Matrix<V> & m){
    this->operator=((*this) & m);
    return *this;
}

template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator&(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.r() != m2.r()) throw std::invalid_argument("Matrices must have same number of rows");

    Matrix<double> ret(m1.r(), m1.c()+m2.c());

    for(uint i=0; i<ret.r(); ++i){
        for(uint j=0; j<ret.c(); ++j){
            if(j < m1.c()) ret(i,j) = m1(i,j);
            else ret(i,j) = m2(i,j-m1.c());
        }
    }

    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator&(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.r() != m2.r()) throw std::invalid_argument("Matrices must have same number of rows");

    Matrix<std::complex<double>> ret(m1.r(), m1.c()+m2.c());

    for(uint i=0; i<ret.r(); ++i){
        for(uint j=0; j<ret.c(); ++j){
            if(j < m1.c()) ret(i,j) = m1(i,j);
            else ret(i,j) = m2(i,j-m1.c());
        }
    }

    return ret;
}


template<typename T>
template <typename U, typename V,
    typename std::enable_if<
        !is_complex<V>::value || 
        (is_complex<U>::value && is_complex<V>::value), 
    int>::type
>
Matrix<T>& Matrix<T>::operator|=(const Matrix<V> & m){
    this->operator=((*this) | m);
    return *this;
}

template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    Matrix<double>
>::type
operator|(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.c() != m2.c()) throw std::invalid_argument("Matrices must have same number of columns");

    Matrix<double> ret(m1.r()+m2.r(), m1.c());

    for(uint i=0; i<ret.r(); ++i){
        for(uint j=0; j<ret.c(); ++j){
            if(i < m1.r()) ret(i,j) = m1(i,j);
            else ret(i,j) = m2(i-m1.r(),j);
        }
    }

    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
operator|(const Matrix<U>& m1, const Matrix<V>& m2){
    if(m1.c() != m2.c()) throw std::invalid_argument("Matrices must have same number of columns");

    Matrix<std::complex<double>> ret(m1.r()+m2.r(), m1.c());

    for(uint i=0; i<ret.r(); ++i){
        for(uint j=0; j<ret.c(); ++j){
            if(i < m1.r()) ret(i,j) = m1(i,j);
            else ret(i,j) = m2(i-m1.r(),j);
        }
    }

    return ret;
}

#pragma endregion concatenate


#pragma region output

template<typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& m){
    os << "Matrix(" << m.r() << "x" << m.c() << ")" << std::endl;
    for(uint i=0; i<m.r(); ++i){
        if(i>0) os << std::endl;
        for(uint j=0; j<m.c(); ++j){
            if(j>0) os << " ";
            os << m(i,j);
        }
    }
    
    return os;
}

std::ostream& operator<<(std::ostream& os, const uu_pair & p){
    os << "Pair(" << p.first << "; " << p.second << ")";
    return os;
}

#pragma endregion output

} // namespace MA

#endif // MA_MATRICES_OP_TPP
