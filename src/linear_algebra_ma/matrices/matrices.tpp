#ifndef MA_MATRICES_TPP
#define MA_MATRICES_TPP

#include "linear_algebra_ma/matrices.hpp"

namespace MA
{

#pragma region static_methods
template<typename T>
T Matrix<T>::rand(){
    return unif(re);
}

template<>
std::complex<double> Matrix<std::complex<double>>::rand(){
    return std::complex<double>(unif(re), unif(re));
}

#pragma endregion static_methods


#pragma region constructor_destructor

template<typename T>
Matrix<T>::Matrix(){
    this->_r = 0;
    this->_c = 0;
    this->_v = nullptr;
}

template<typename T>
Matrix<T>::Matrix(const uint & r, const uint & c) {
    this->_r = r;
    this->_c = c;
    this->_v = new T[this->size()]();
}

template<typename T>
Matrix<T>::Matrix(const uint & r, const uint & c, std::vector<T> v){
    this->_r = r;
    this->_c = c;
    this->_v = new T[this->size()];
    this->set(v);
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
Matrix<T>::Matrix(const uint & r, const uint & c, std::vector<U> v) {
    this->_r = r;
    this->_c = c;
    this->_v = new T[this->size()];
    this->set(v);   
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T> & m){
    this->_r = m.r();
    this->_c = m.c();
    this->_v = new T[m.size()];
    std::copy(m.v(), m.v() + m.size(), this->_v);
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
Matrix<T>::Matrix(const Matrix<U> & m){
    this->_r = m.r();
    this->_c = m.c();
    this->_v = new T[m.size()];
    for(uint i=0; i<m.size(); ++i){
        this->_v[i] = m.v()[i];
    }
}

template<typename T>
Matrix<T>::Matrix(Matrix<T> && m) noexcept
: Matrix()
{
    swap(*this, m);
}

template<typename T>
Matrix<T>::~Matrix() {
    if(this->_v != nullptr) delete[] this->_v;
}

#pragma endregion constructor_destructor


#pragma region getter

template<typename T>
uint Matrix<T>::r() const{return _r;}
template<typename T>
uint Matrix<T>::c() const{return _c;}
template<typename T>
uu_pair Matrix<T>::shape() const{return uu_pair{_r, _c};}
template<typename T>
uint Matrix<T>::size() const{return _r * _c;}
template<typename T>
T const * Matrix<T>::v() const{return _v;}

template<typename T>
Matrix<T> Matrix<T>::diag() const{
    uint dim = std::min(_c, _r);
    Matrix v = Matrix(dim,1);

    for(uint i=0; i<dim; ++i){
        v(i) = this->operator()(i,i);
    }

    return v;
}

template<typename T>
template <typename U, 
    typename std::enable_if<is_complex<U>::value, int>::type
>
Matrix<double> Matrix<T>::real() const{
    std::vector<double> v;
    for(uint i=0; i<this->size(); ++i){
        v.push_back(this->_v[i].real());
    }
    return Matrix<double>(_r, _c, v);
}

template<typename T>
template <typename U, 
    typename std::enable_if<is_complex<U>::value, int>::type
>
Matrix<double> Matrix<T>::imag() const{
    std::vector<double> v;
    for(uint i=0; i<this->size(); ++i){
        v.push_back(this->_v[i].imag());
    }
    return Matrix<double>(_r, _c, v);
}

#pragma endregion getter


#pragma region setter

template<typename T>
void Matrix<T>::set(std::vector<T> v){
    if(v.size() < this->size()) throw std::out_of_range("Not enough values in v to init the matrix");

    std::copy(v.begin(), v.begin() + this->size(), this->_v);
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
void Matrix<T>::set(std::vector<U> v){
    if(v.size() < this->size()) throw std::out_of_range("Not enough values in v to init the matrix");

    for(uint i=0; i<this->size(); ++i){
        this->_v[i] = v[i];
    }
}

template<typename T>
void Matrix<T>::set(uu_pair rs, uu_pair cs, std::vector<T> v){
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(v.size() < nRows*nCols) throw std::out_of_range("Given vector doesn't have enough elements");

    for (uint i = 0; i < nRows; ++i){
        std::copy(v.begin() + i * nCols, v.begin() + (i+1) * nCols, this->_v + (rs.first + i) * _c + cs.first);
    }
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
void Matrix<T>::set(uu_pair rs, uu_pair cs, std::vector<U> v){
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(v.size() < nRows*nCols) throw std::out_of_range("Given vector doesn't have enough elements");

    for (uint i = 0; i < nRows; ++i){
        for (uint j = 0; j < nCols; ++j){
            this->operator()(i + rs.first, j + cs.first) = v[i * nCols + j];
        }
    }
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
void Matrix<T>::set(uint r, uint c, U x){
    this->operator()(r,c) = x;
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
void Matrix<T>::set(uu_pair rs, uint c, Matrix<U> m){
    if(rs.second >= this->_r) throw std::out_of_range("Row index out of range");
    if(c >= this->_c) throw std::out_of_range("Col index out of range");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second");

    // allow to use {} for full row selection
    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;

    if(m.shape() != uu_pair{rs.second - rs.first + 1, 1}) throw std::invalid_argument("Given matrix's shape does not match");

    for(uint i=0; i<m._r; ++i){
        this->operator()(i + rs.first, c) = m(i,0);
    }
}

template<typename T>
void Matrix<T>::set(uint r, uu_pair cs, Matrix<T> m){
    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw std::out_of_range("Col index out of range");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    // allow to use {} for full col selection
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    if(m.shape() != uu_pair{1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");

    for(uint j=0; j<m._c; ++j){
        std::copy(m._v, m._v + m.size(), this->_v + _c*r + cs.first);
    }
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
void Matrix<T>::set(uint r, uu_pair cs, Matrix<U> m){
    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw std::out_of_range("Col index out of range");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    // allow to use {} for full col selection
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    if(m.shape() != uu_pair{1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");

    for(uint j=0; j<m.c(); ++j){
        this->operator()(r, j + cs.first) = m(r,j);
    }
}

template<typename T>
void Matrix<T>::set(uu_pair rs, uu_pair cs, Matrix<T> m){
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    // allow to use {} for full row/col selection
    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    if(m.shape() != uu_pair{rs.second - rs.first + 1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");;

    for (uint i = 0; i < m.r(); ++i){
        std::copy(m._v + i * m.c(), m._v + (i+1) * m.c(), this->_v + (rs.first + i) * _c + cs.first);
    }
}

template<typename T>
template<typename U,
    typename std::enable_if<
        !is_complex<U>::value ||
        (is_complex<T>::value && is_complex<U>::value), int
    >::type
>
void Matrix<T>::set(uu_pair rs, uu_pair cs, Matrix<U> m){
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    // allow to use {} for full row/col selection
    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    if(m.shape() != uu_pair{rs.second - rs.first + 1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");;

    for (uint i = 0; i < m.r(); ++i){
        for (uint j = 0; j < m.c(); ++j){
            this->operator()(i + rs.first, j + cs.first) = m(i,j);
        }
    }
}

template<typename T>
Matrix<T> Matrix<T>::reshape(const uint & r, const uint & c) const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    if(r*c != this->size()) throw std::invalid_argument("New matrix size must match the current one");

    Matrix ret(r,c);
    std::copy(_v, _v + this->size(), ret._v);
    return ret;
}

#pragma endregion setter


#pragma region swap

template<typename T>
void swap(Matrix<T> & m1, Matrix<T> & m2){
    using std::swap;
    // swap attributes
    std::swap(m1._r, m2._r);
    std::swap(m1._c, m2._c);
    std::swap(m1._v, m2._v);
}

template<typename T>
void Matrix<T>::swap_rows(const uint & r1, const uint & r2){
    if(r1 >= _r || r2 >= _r) throw std::out_of_range("Given parameters exceed the matrix rows indeces");

    // swap not necessary
    if(r1 == r2) return;

    // extract r1
    Matrix tmp = this->operator()(r1, ALL);
    // substitute r2 into r1
    this->set(r1, ALL, this->operator()(r2, ALL));
    // substitute r1 into r2
    this->set(r2, ALL, tmp);
}

template<typename T>
void Matrix<T>::swap_cols(const uint & c1, const uint & c2){
    if(c1 >= _c || c2 >= _c) throw std::out_of_range("Given parameters exceed the matrix columns indeces");

    // swap not necessary
    if(c1 == c2) return;

    // extract c1
    Matrix tmp = this->operator()(ALL, c1);
    // substitute c2 into c1
    this->set(ALL, c1, this->operator()(ALL, c2));
    // substitute c1 into c2
    this->set(ALL, c2, tmp);
}

#pragma endregion swap


#pragma region vector

template<typename T>
Matrix<T> Matrix<T>::to_c_vec() const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    return this->reshape(this->size(), 1);
}

template<typename T>
Matrix<T> Matrix<T>::to_r_vec() const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    return this->reshape(1, this->size());
}

template<typename U, typename V>
typename std::enable_if<
    !is_complex<U>::value && !is_complex<V>::value, 
    double
>::type
dot(const Matrix<U> & v1, const Matrix<V> & v2){
    if(! v1.is_vec()) throw std::invalid_argument("Obj 1 is not a vector");
    if(! v2.is_vec()) throw std::invalid_argument("Obj 2 is not a vector");
    if(v1.size() != v2.size()) throw std::invalid_argument("Vectors length don't match");

    double ret = 0;
    for(uint i=0; i<v1.size(); ++i) ret += v1(i) * v2(i);
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    std::complex<double>
>::type
dot(const Matrix<U> & v1, const Matrix<V> & v2){
    if(! v1.is_vec()) throw std::invalid_argument("Obj 1 is not a vector");
    if(! v2.is_vec()) throw std::invalid_argument("Obj 2 is not a vector");
    if(v1.size() != v2.size()) throw std::invalid_argument("Vectors length don't match");

    std::complex<double> ret = 0.0;
    for(uint i=0; i<v1.size(); ++i) ret += v1(i) * v2(i);
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<double>
>::type
cross(const Matrix<U> & v1, const Matrix<V> & v2){
    if(! v1.is_vec()) throw std::invalid_argument("Obj 1 is not a vector");
    if(! v2.is_vec()) throw std::invalid_argument("Obj 2 is not a vector");
    if(v1.size() != 3 || v1.size() != v2.size()) throw std::invalid_argument("Cross product defined only for 3d vectors");

    Matrix<double> ret(3,1);
    ret(0) =   (v1(1) * v2(2)) - (v1(2) * v2(1));
    ret(1) = -((v1(0) * v2(2)) - (v1(2) * v2(0)));
    ret(2) =   (v1(0) * v2(1)) - (v1(1) * v2(0));
    return ret;
}

template<typename U, typename V>
typename std::enable_if<
    is_complex<U>::value || is_complex<V>::value, 
    Matrix<std::complex<double>>
>::type
cross(const Matrix<U> & v1, const Matrix<V> & v2){
    if(! v1.is_vec()) throw std::invalid_argument("This obj is not a vector");
    if(! v2.is_vec()) throw std::invalid_argument("Given obj is not a vector");
    if(v1.size() != 3 || v1.size() != v2.size()) throw std::invalid_argument("Cross product defined only for 3d vectors");

    Matrix<std::complex<double>> ret(3,1);
    ret(0) =   (v1(1) * v2(2)) - (v1(2) * v2(1));
    ret(1) = -((v1(0) * v2(2)) - (v1(2) * v2(0)));
    ret(2) =   (v1(0) * v2(1)) - (v1(1) * v2(0));
    return ret;
}

template<typename T>
double Matrix<T>::norm2() const{
    if(this->_r == 1 || this->_c == 1){
        double ret = 0;
        for(uint i=0; i< this->_r+this->_c-1; ++i){
            double tmp = abs(_v[i]);
            ret += tmp * tmp;
        }
        return sqrt(ret);
    }
    else{
        throw std::invalid_argument("Norm only appliable to row/column vectors");
    }
}

template<typename T>
Matrix<T> Matrix<T>::normalize() const{
    return Matrix<T>(*this) / this->norm2();
}

template<typename T>
Matrix<T>& Matrix<T>::normalize_self(){
    return this->operator/=(this->norm2());
}

#pragma endregion vector


#pragma region checks

template<typename T>
bool Matrix<T>::is_vec() const{
    return this->_r == 1 || this->_c == 1;
}

template<typename T>
bool Matrix<T>::is_sing() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");
    return this->det() == 0.0;
}

template<typename T>
bool Matrix<T>::is_upper_triang() const{
    uint n = std::min(_r, _c);
    for(uint i=1; i<n; ++i){
        for(uint j=0; j<i; ++j){
            if(abs(this->operator()(i,j)) > Matrix<T>::epsilon) return false;
        }
    }

    return true;
}
    
template<typename T>
bool Matrix<T>::is_lower_triang() const{
    uint n = std::min(_r, _c);
    for(uint j=1; j<n; ++j){
        for(uint i=0; i<j; ++i){
            if(abs(this->operator()(i,j)) > Matrix<T>::epsilon) return false;
        }
    }

    return true;
}

template<typename T>
bool Matrix<T>::is_upper_hessenberg() const{
    uint n = std::min(_r, _c);
    for(uint i=2; i<n; ++i){
        for(uint j=0; j<i-1; ++j){
            if(abs(this->operator()(i,j)) > Matrix<T>::epsilon) return false;
        }
    }

    return true;
}

template<typename T>
bool Matrix<T>::is_lower_hessenberg() const{
    uint n = std::min(_r, _c);
    for(uint j=2; j<n; ++j){
        for(uint i=0; i<j-1; ++i){
            if(abs(this->operator()(i,j)) > Matrix<T>::epsilon) return false;
        }
    }

    return true;
}

#pragma endregion checks


#pragma region matrix_operations

template<typename T>
Matrix<T> Matrix<T>::t() const{
    Matrix<T> ret(this->_c, this->_r);
    for (uint i = 0; i<this->_r; ++i) {
        for (uint j = 0; j<this->_c; ++j) {
            ret._v[j * ret._c + i] = this->_v[i * this->_c + j];
        }
    }
    return ret;
}

template <>
Matrix<std::complex<double>> Matrix<std::complex<double>>::t() const{
    Matrix<std::complex<double>> ret(this->_c, this->_r);
    for (uint i = 0; i<this->_r; ++i) {
        for (uint j = 0; j<this->_c; ++j) {
            ret._v[j * ret._c + i] = std::conj(this->_v[i * this->_c + j]);
        }
    }
    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::submat_del(const uint & p, const uint & q) const{
    if(this->_r <= p) throw std::invalid_argument("p greater than matrix rows");
    if(this->_c <= q) throw std::invalid_argument("q greater than matrix cols");

    uint i = 0, j = 0;

    Matrix<T> ret = Matrix(this->_r-1,this->_c-1);

    // Looping for each element of the matrix
    for (uint row = 0; row < this->_r; row++){

        for (uint col = 0; col < this->_c; col++){

            // Copying into result matrix only those element
            // which are not in given row and column
            if (row != p && col != q){
                ret(i,j) = this->operator()(row,col);
                ++j;
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == ret._c){
                    j = 0;
                    ++i;
                }
            }
        }
    }

    return ret;
}

template<typename T>
T Matrix<T>::det() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    // empty matrix: return 0
    if (this->_r == 0){
        return 0;
    } 
    // matrix contains single element
    if (this->_r == 1){
        return this->_v[0];
    } 
    // matrix is 2x2
    else if (this->_r == 2){
        return this->_v[0]*this->_v[3] - this->_v[1]*this->_v[2];
    }
    // matrix is bigger then 2x2
    else{
        Matrix<T> L, U; 
        Matrix<double> P;
        uint n_swaps = this->lup_dec(L,U,P);
        T determinant = n_swaps % 2 == 0 ? 1.0 : -1.0;
        for(uint i=0; i<U.c(); ++i) determinant *= U(i,i);
        return determinant;
    }
}


template<typename T>
T Matrix<T>::minor(const uint & p, const uint & q) const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    return this->submat_del(p,q).det();
}

template<typename T>
T Matrix<T>::cof(const uint & p, const uint & q) const{
    return (((p+q) % 2) == 0 ? 1.0 : -1.0) * this->minor(p, q);
}

template<typename T>
Matrix<T> Matrix<T>::cof_mat() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    Matrix<T> ret(this->_r, this->_c);
    for(uint i=0; i<this->_r; ++i) for(uint j=0; j<this->_c; ++j){
        ret(i,j) = this->cof(i,j);
    }

    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::adj() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    // same formula of cofactor matrix, but inverting indeces to get the transpose
    Matrix<T> ret(this->_r, this->_c);
    for(uint i=0; i<this->_r; ++i) for(uint j=0; j<this->_c; ++j){
        ret(j,i) = this->cof(i,j);
    }

    return ret;
}

// Matrix Matrix<T>::inv() const{
//     if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

//     return matrix_l_divide(*this, IdMat(_r));
// }

// Matrix Matrix<T>::pinv_left() const{
//     if(_c > _r) throw std::invalid_argument("Matrix must be full column rank");
//     try{
//         Matrix tmp = this->t() * *this;
//         double det = tmp.det();
//         if(abs(det) < Matrix<T>::epsilon) {
//             std::cout << "Matrix is bad conditioned, det = " << det << "; this may result in bad numerical result" << std::endl;
//         }
//         return tmp.inv() * this->t();
//     }
//     catch(const std::runtime_error& rte){
//         std::cerr << rte.what() << '\n';
//         throw std::runtime_error("Matrix (m.t * m) not invertible");
//     }
//     catch(const std::exception& e){
//         std::cerr << e.what() << '\n';
//         throw std::runtime_error("Unknown error occured");
//     }
// }

// Matrix Matrix<T>::pinv_right() const{
//     if(_r > _c) throw std::invalid_argument("Matrix must be full row rank");
//     try{
//         Matrix tmp = (*this * this->t());
//         double det = tmp.det();
//         if(abs(det) < Matrix<T>::epsilon) {
//             std::cout << "Matrix is bad conditioned, det = " << det << "; this may result in bad numerical result" << std::endl;
//         }
//         return this->t() * tmp.inv();
//     }
//     catch(const std::runtime_error& rte){
//         std::cerr << rte.what() << '\n';
//         throw std::runtime_error("Matrix (m.t * m) not invertible");
//     }
//     catch(const std::exception& e){
//         std::cerr << e.what() << '\n';
//         throw std::runtime_error("Unknown error occured");
//     }
// }

#pragma endregion matrix_operations


#pragma region decomposition_methods

template<typename T>

Matrix<T> Matrix<T>::householder_v() const{
    if(!this->is_vec()) throw std::invalid_argument("The object must be a vector");

    using namespace std::complex_literals;
    Matrix<T> u(*this);
    T alpha;

    if constexpr (is_complex<T>::value) alpha = -exp(arg(u(0))*1i) * u.norm2();
    else alpha = -copysign(u.norm2(), u(0));

    u(0) -= alpha;
    return u.normalize();
}

template<typename T>
void Matrix<T>::qr_dec(Matrix<T> & Q, Matrix<T> & R) const{
    // init matrices
    uint n = std::min<double>(_r, _c);
    Q = IdMat(_r);
    R = *this;
    
    for(uint i=0; i<n-1; ++i){
        //compute vk
        Matrix<T> v = R({i, _r-1}, i).householder_v();
        Matrix<T> v_t = v.t();

        // Update R: HR = R - 2v * (v' * R)
        for(uint j=0; j<_c; ++j){
            T tmp = 0.0;
            for(uint k=0; k<v.size(); ++k) tmp += v_t(k) * R(i+k,j);
            for(uint k=0; k<v.size(); ++k) R(i+k,j) -= 2.0 * v(k) * tmp;
        }
        // Update Q: QH = Q - (Q * v) * 2v'
        for(uint j=0; j<_r; ++j){
            T tmp = 0.0;
            for(uint k=0; k<v.size(); ++k) tmp += v(k) * Q(j,i+k);
            for(uint k=0; k<v.size(); ++k) Q(j,i+k) -= 2.0 * v_t(k) * tmp;
        }
    }
}

template<typename T>
void Matrix<T>::qrp_dec(Matrix<T> & Q, Matrix<T> & R, Matrix<double> & P) const{

    uint n = std::min<double>(_r, _c);

    // init matrices
    Q = IdMat(_r);
    R = *this;
    P = IdMat(_c);
    
    for(uint i=0; i<n-1; ++i){ // main loop

        // find column with largest norm
        uint index = 0;
        double max_norm = -1;
        for(uint k=i; k<n; ++k){
            double norm = R({i, _r-1}, k).norm2();
            if(norm > max_norm){
                max_norm = norm;
                index = k;
            }
        }

        // swap columns
        R.swap_cols(i, index);
        P.swap_cols(i, index);

        //compute vk
        Matrix<T> v = R({i, _r-1}, i).householder_v();
        Matrix<T> v_t = v.t();

        // Update R: HR = R - 2v * (v' * R)
        for(uint j=0; j<_c; ++j){
            T tmp = 0.0;
            for(uint k=0; k<v.size(); ++k) tmp += v_t(k) * R(i+k,j);
            for(uint k=0; k<v.size(); ++k) R(i+k,j) -= 2.0 * v(k) * tmp;
        }
        // Update Q: QH = Q - (Q * v) * 2v'
        for(uint j=0; j<_r; ++j){
            T tmp = 0.0;
            for(uint k=0; k<v.size(); ++k) tmp += v(k) * Q(j,i+k);
            for(uint k=0; k<v.size(); ++k) Q(j,i+k) -= 2.0 * v_t(k) * tmp;
        }
    }
}

template<typename T>
uint Matrix<T>::lup_dec(Matrix<T> & L, Matrix<T> & U, Matrix<double> & P) const{
    if(_r != _c) throw std::invalid_argument("The matrix must be square");

    L = IdMat(_r);
    U = *this;
    P = IdMat(_r);
    uint ret = 0;

    for (uint i = 0; i < _c; ++i) { // scroll columns
        // pivoting
        double u_max = 0.0;
        uint max_index = i;
        for (uint j = i; j < _r; ++j){ // scroll rows
            // find max value
            if(u_max < abs(U(j,i))){
                max_index = j;
                u_max = abs(U(j,i));
            }
        }
        // if max value is zero we can skip this iteration
        if(u_max == 0.0) continue;

        // swap rows
        if(max_index != i){
            P.swap_rows(i, max_index);
            U.swap_rows(i, max_index);
            ret++;
        }

        // elimination
        for (uint j = i+1; j < _r; ++j){ // scroll rows
            // compute l(j,i)
            U(j,i) = U(j,i) / U(i,i);
            // apply elimination over all elements of the row
            for (uint k = i+1; k < _c; ++k){ // scroll columns
                U(j,k) -= U(j,i) * U(i,k);
            }
        }    
    }

    // compose matrices
    for(uint i=1; i<_r; ++i) { //scroll rows
        for(uint j=0; j<i; ++j){ //scroll columns
            L(i,j) = U(i,j);
            U(i,j) = 0.0;
        }
    }

    return ret;
}

template<typename T>
void Matrix<T>::hessenberg_dec(Matrix<T> & Q, Matrix<T> & H) const{
    if(_c != _r) throw std::invalid_argument("Matrix must be square");

    H = *this;
    Q = IdMat(_r);
    
    using namespace std;

    for (uint i=0; i<_r-1; ++i){
        // compute vk
        Matrix<T> v = H({i+1, _r-1}, i).householder_v();

        // compute U matrix
        Matrix<T> U = IdMat(_r);
        U.set({i+1, _r-1},{i+1, _r-1}, IdMat(v.size()) - 2 * v * v.t());

        // // Update R: HR = R - 2v * (v' * R)
        // for(uint j=0; j<_c; ++j){
        //     T tmp = 0.0;
        //     for(uint k=0; k<v.size(); ++k) tmp += v_t(k) * R(i+k,j);
        //     for(uint k=0; k<v.size(); ++k) R(i+k,j) -= 2.0 * v(k) * tmp;
        // }
        // // Update Q: QH = Q - (Q * v) * 2v'
        // for(uint j=0; j<_r; ++j){
        //     T tmp = 0.0;
        //     for(uint k=0; k<v.size(); ++k) tmp += v(k) * Q(j,i+k);
        //     for(uint k=0; k<v.size(); ++k) Q(j,i+k) -= 2.0 * v_t(k) * tmp;
        // }

        // // Update H: UH = H - 2v * (v' * H) -> H'
        // for(uint j=0; j<_c; ++j){
        //     T tmp = 0.0;
        //     for(uint k=0; k<v.size(); ++k) tmp += v_t(k) * H(i+k,j);
        //     for(uint k=0; k<v.size(); ++k) H(i+k,j) -= 2.0 * v(k) * tmp;
        // }
        // // Update H': H'U.t = H - (H * v) * 2v'

        H = U * H * U;
        Q = Q * U;
    }
}

#pragma endregion decomposition_methods


#pragma region ls_solution

// Matrix Matrix<T>::backward_sub(Matrix const & U, Matrix const & B){
//     if(U.c() != U.r()) throw std::invalid_argument("Coefficient matrix U must be square");
//     if(B.r() != U.c()) throw std::invalid_argument("Rows of B must be equal to the columns of U");
//     // check all U diagonal elements are different from zero
//     for(uint i=0; i<U.c(); ++i) if(U(i,i) == 0){
//         throw std::runtime_error("System of equation is underdetermined");
//     }

//     Matrix res(U.r(), B.c());

//     double tmp;
//     // printf("start i loop\n");
//     for(uint i=0; i<B.c(); ++i){ //col of res == col of B
//         // printf("start j loop with i:%d\n", i);
//         for(int j=U.r()-1; j>=0; j--){ //row of res == row of U == row of B
//             tmp = 0;
//             // printf("start k loop with i:%d, j:%d\n", i, j);
//             for(int k=U.r()-1; k>j; k--){ //col of U = row of res
//                 tmp += U(j,k) * res(k,i);
//                 // printf("i:%d j:%d k:%d\n",i,j,k);
//             }
            
//             res(j, i) = (B(j,i) - tmp) / U(j,j);
//             // printf("%f\t%f\n",tmp,res[j*B.c() + i]);
//         }
//         // printf("\n");
//     }
    
//     return res;
// }

// Matrix Matrix<T>::forward_sub(Matrix const & L, Matrix const & B){
//     if(L.c() != L.r()) throw std::invalid_argument("Coefficient matrix L must be square");
//     if(B.r() != L.c()) throw std::invalid_argument("Rows of B must be equal to the columns of L");
//     // check all L diagonal elements are different from zero
//     for(uint i=0; i<L.c(); ++i) if(L(i,i) == 0){
//         throw std::runtime_error("System of equation is underdetermined");
//     }

//     Matrix res(L.r(), B.c());
    
//     double tmp;
//     for(uint j=0; j<L.r(); ++j){ //row of res == row of L == row of B
//         // printf("j: %d\n", j);
//         for(uint i=0; i<B.c(); ++i){ //col of res == col of B
//             // printf("i: %d\n", i);
//             tmp = 0;
//             for(uint k=0; k<j; ++k){ //col of L = row of res
//                 // printf("k: %d\n", k);
//                 tmp += L(j,k) * res(k,i);
//             }
//             // printf("%f\n", tmp);
//             res(j,i) = (B(j,i) - tmp) / L(j,j);
//         }
//     }
    
//     return res;
// }

// Matrix Matrix<T>::matrix_l_divide(Matrix const & A, Matrix const & B){
//     if(A.r() != B.r()){
//         throw std::invalid_argument("Rows of B (currently " + std::to_string(B.r()) + ")" +
//         "must be equal the rows of A (currently " + std::to_string(A.r()) + ")");
//     }

//     if(A.c() == A.r()){ //square A
//         // printf("matrix A is square, using LU decomposition\n");
//         Matrix L, U, P, res_tmp;

//         A.lup_dec(L, U, P);
//         res_tmp = forward_sub(L, P*B);
//         return backward_sub(U, res_tmp);
//     }
//     else{
//         // printf("matrix A is not square, using QR decomposition\n");
//         if(A.r() < A.c()) throw std::invalid_argument("System is underdetermined");
//         Matrix Q,R,P, res;

//         A.qrp_dec(Q,R,P);
//         res = P * backward_sub(R({0, R.c()-1}, ALL), (Q.t() * B)({0, R.c()-1}, ALL));
//         return res;
//     }

// }

// Matrix Matrix<T>::solve_ls(Matrix const & A, Matrix const & B){
//     return Matrix<T>::matrix_l_divide(A,B);
// }

// Matrix Matrix<T>::matrix_r_divide(Matrix const & B, Matrix const & A){
//     if(A.c() != B.c()){
//         throw std::invalid_argument("Columns of B (currently " + std::to_string(B.c()) + ") " +
//         "must be equal the Columns of A (currently " + std::to_string(A.c()) + ")");
//     }
//     return matrix_l_divide(A.t(), B.t()).t();
// }

#pragma endregion ls_solution


#pragma region eigen

// void Matrix<T>::eigen_QR(Matrix & D, Matrix & V, uint max_iterations, double tolerance) const{
//     if(_c != _r) throw std::invalid_argument("Matrix must be square");

//     Matrix Q, R, P;
//     uint k;
//     double sum;
//     this->hessenberg_dec(Q,D);

//     for(k=0; k<max_iterations; ++k){
//         Matrix tmp = D.diag();
//         // D.qrp_dec(Q,R,P);
//         // D = R*Q;
//         D.qr_dec(Q,R);
//         D = R*Q;

//         // check for convergence
//         // double max_variation = -1;
//         // tmp = tmp - D.diag();
//         // for(uint i=0; i<_r; ++i){
//         //     max_variation = std::max(max_variation, abs(tmp(i) / D(i)));
//         // }
//         sum = 0;
//         for(uint i=1; i<_r; ++i){
//             for(uint j=0; j<i; ++j){
//                 sum += abs(D(i,j));
//             }
//         }
//         if(sum < tolerance) break;
//     }

//     // std::cout << "QR algorithm steps: " << k << ", residuals: " << sum << std::endl;
//     // std::cout << D << std::endl;
//     D = D.diag();

//     V = Matrix(_r, _r);
//     for(uint i=0; i<_r; ++i){
//         Matrix M = *this - D(i) * IdMat(_r);
//         Matrix eigen_vec = Matrix(1,1,{1}) | solve_ls(M(ALL, {1,_r-1}), - M(ALL, 0));
//         V.set(ALL, i, eigen_vec.normalize());
//     }

// }

template<typename T>
Matrix<T> Matrix<T>::implicit_double_QR_step() const{
    // for easier notation
    uint n =  this->_r;
    Matrix A = *this;

    Matrix y,u,v;
    double tau, tmp;

    // -- compute and apply Q1
    // compute first column of B
    y = Matrix(3,1,{
        (   (A(0,0) - A(n-2,n-2)) * 
            (A(0,0) - A(n-1, n-1)) - 
            A(n-1,n-2) * A(n-2,n-1)
        ) / A(1,0) + A(0,1),
        A(0,0) + A(1,1) - 
            A(n-2,n-2) - A(n-1,n-1),
        A(2,1)
    });
    // compute reflector
    tau = copysign(y.norm2(), y(0));
    u = Matrix(3,1,{1, y(1) / (y(0) + tau), y(2) / (y(0) + tau) });
    v = (y(0) / tau + 1) * u;
    // apply reflection B -> QB = B - v * (u.t * B)
    for(uint j=0; j<n; j++){
        tmp = 0;
        for(uint k=0; k<3; k++) tmp += u(k) * A(k, j);
        for(uint k=0; k<3; k++) A(k, j) -= tmp * v(k);
    }
    // apply reflection C -> CQ = C - (C * u) * v.t
    for(uint j=0; j<n; j++){
        tmp = 0;
        for(uint k=0; k<3; k++) tmp += v(k) * A(j, k);
        for(uint k=0; k<3; k++) A(j, k) -= tmp * u(k);
    }
    // std::cout << "after applying Q1: " << A << std::endl << std::endl;

    // -- move the bulge
    for(uint i=1; i<n-2; ++i){
        // std::cout << "acting on column # " << i-1 << std::endl;
        // extract column to transform in  [x 0 0]
        y = A({i, i+2}, i-1);
        // compute reflector
        tau = copysign(y.norm2(), y(0));
        u = Matrix(3,1,{1, y(1) / (y(0) + tau), y(2) / (y(0) + tau)});
        v = (y(0) / tau + 1) * u;
        // apply reflection B -> QB = B - v * (u.t * B)
        for(uint j=0; j<n; j++){
            tmp = 0;
            for(uint k=0; k<3; k++) tmp += u(k) * A(i+k, j);
            for(uint k=0; k<3; k++) A(i+k, j) -= tmp * v(k);
        }
        // apply reflection C -> CQ = C - (C * u) * v.t
        for(uint j=0; j<n; j++){
            tmp = 0;
            for(uint k=0; k<3; k++) tmp += v(k) * A(j, i+k);
            for(uint k=0; k<3; k++) A(j, i+k) -= tmp * u(k);
        }
        // std::cout << "moving the bulge: " << A << std::endl << std::endl;

    }

    // std::cout << "last iteration" << std::endl;
    // -- implement step for n-1 with only 2 last rows
    y = A({n-2, n-1}, n-3);
    tau = copysign(y.norm2(), y(0));
    u = Matrix(2,1,{1, y(1) / (y(0) + tau)});
    v = (y(0) / tau + 1) * u;
    // apply reflection B -> QB = B - v * (u.t * B)
    for(uint j=0; j<n; j++){
        tmp = 0;
        for(uint k=0; k<2; k++) tmp += u(k) * A(n-2+k, j);
        for(uint k=0; k<2; k++) A(n-2+k, j) -= tmp * v(k);
    }
    // apply reflection C -> CQ = C - (C * u) * v.t
    for(uint j=0; j<n; j++){
        tmp = 0;
        for(uint k=0; k<2; k++) tmp += v(k) * A(j, n-2+k);
        for(uint k=0; k<2; k++) A(j, n-2+k) -= tmp * u(k);
    }
    // std::cout << "end of QR step: " << A << std::endl << std::endl;

    return A;
}

// void Matrix<T>::eigen_implicit_QR(Matrix & D, Matrix & V, uint max_iterations, double tolerance) const{
//     if(_c != _r) throw std::invalid_argument("Matrix must be square");

//     Matrix Q, R, P;
//     // obtain an hessenberg matrix
//     this->hessenberg_dec(Q,D);

//     // init shift
//     double mu = D(_r-1);
//     uint k;
//     double sum;

//     for(k=0; k<max_iterations; ++k){
//         // implicit double QR step
        


//         Matrix tmp = D.diag();
        
//         // apply shift
//         (D - mu * IdMat(_r)).qr_dec(Q,R);
//         // shift backward
//         D = R*Q + mu * IdMat(_r);
//         // update mu
//         Matrix v = Q(ALL, _c-1);
//         mu = (v.t() * D * v).operator()(0) / (v.t() * v).operator()(0);

//         // check for convergence
//         // double max_variation = -1;
//         // tmp = tmp - D.diag();
//         // for(uint i=0; i<_r; ++i){
//         //     max_variation = std::max(max_variation, abs(tmp(i) / D(i)));
//         // }
//         sum = 0;
//         for(uint i=1; i<_r; ++i){
//             for(uint j=0; j<i; ++j){
//                 sum += abs(D(i,j));
//             }
//         }
//         if(sum < tolerance) break;
//     }

//     std::cout << "QR with shift algorithm steps: " << k << ", residuals: " << sum << std::endl;
//     std::cout << D << std::endl;
//     D = D.diag();

//     V = Matrix(_r, _r);
//     for(uint i=0; i<_r; ++i){
//         Matrix M = *this - D(i) * IdMat(_r);
//         Matrix eigen_vec = Matrix(1,1,{1}) | solve_ls(M(ALL, {1,_r-1}), - M(ALL, 0));
//         V.set(ALL, i, eigen_vec.normalize());
//     }

// }

#pragma endregion eigen


#pragma region special_constructors

Matrix<double> IdMat(const uint & dim){
    Matrix ret = Matrix(dim,dim);
    for(uint i=0; i<dim; ++i){
        ret(i,i) = 1;
    }
    return ret;
}

// Matrix<double> Ones(const uint & r, const uint & c){
//     return Matrix(r,c)+1;
// }

template<typename T>
Matrix<T> RandMat(const uint & r, const uint & c){
    Matrix<T> ret(r,c);
    for(uint i=0; i<r; ++i){
        for(uint j=0; j<c; ++j){
            ret(i,j) = Matrix<T>::rand();
        }
    }

    return ret;
}

template<typename T>
Matrix<T> diag(std::vector<T> v){
    Matrix ret = Matrix(v.size(), v.size());
    
    for(uint i=0; i<v.size(); ++i){
        ret(i,i) = v[i];
    }

    return ret;
}

template<typename T>
Matrix<T> diag(const Matrix<T>& v){
    if(! v.is_vec()) throw std::invalid_argument("Input matrix must be a vector");

    Matrix<T> ret(v.size(), v.size());
    
    for(uint i=0; i<v.size(); ++i){
        ret(i,i) = v(i);
    }

    return ret;
}

#pragma endregion special_constructors

} // namespace MA

#endif // MA_MATRICES_TPP
