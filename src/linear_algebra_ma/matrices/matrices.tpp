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
Matrix<T>::Matrix(uint r, uint c) {
    this->_r = r;
    this->_c = c;
    this->_v = new T[this->size()]();
}

template<typename T>
Matrix<T>::Matrix(uint r, uint c, std::vector<T> v){
    this->_r = r;
    this->_c = c;
    this->_v = new T[this->size()];
    this->set(v);
}

template<typename T>
template<typename U, typename>
Matrix<T>::Matrix(uint r, uint c, std::vector<U> v) {
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
template<typename U, typename>
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
T& Matrix<T>::at(uint r, uint c){
    return this->operator()(r,c);
}

template<typename T>
const T& Matrix<T>::at(uint r, uint c) const{
    return this->operator()(r,c);
}

template<typename T>
T& Matrix<T>::at(uint i){
    return this->operator()(i);
}

template<typename T>
const T& Matrix<T>::at(uint i) const{
    return this->operator()(i);
}

template<typename T>
Matrix<T> Matrix<T>::at(uu_pair rs, uint c) const{
    return this->operator()(rs,c);
}

template<typename T>
Matrix<T> Matrix<T>::at(uint r, uu_pair cs) const{
    return this->operator()(r,cs);
}

template<typename T>
Matrix<T> Matrix<T>::at(uu_pair rs, uu_pair cs) const{
    return this->operator()(rs,cs);
}

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
        v(i) = this->at(i,i);
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
Matrix<T>& Matrix<T>::set(std::vector<T> v){
    if(v.size() < this->size()) throw std::out_of_range("Not enough values in v to init the matrix");

    std::copy(v.begin(), v.begin() + this->size(), this->_v);
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::set(std::vector<U> v){
    if(v.size() < this->size()) throw std::out_of_range("Not enough values in v to init the matrix");

    for(uint i=0; i<this->size(); ++i){
        this->_v[i] = v[i];
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::set(uu_pair rs, uu_pair cs, std::vector<T> v){
    // allow to use -1 to indicate end of row/column
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;

    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(v.size() < nRows*nCols) throw std::out_of_range("Given vector doesn't have enough elements");

    for (uint i = 0; i < nRows; ++i){
        std::copy(v.begin() + i * nCols, v.begin() + (i+1) * nCols, this->_v + (rs.first + i) * _c + cs.first);
    }
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::set(uu_pair rs, uu_pair cs, std::vector<U> v){
    // allow to use -1 to indicate end of row/column
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;

    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(v.size() < nRows*nCols) throw std::out_of_range("Given vector doesn't have enough elements");

    for (uint i = 0; i < nRows; ++i){
        for (uint j = 0; j < nCols; ++j){
            this->at(i + rs.first, j + cs.first) = v[i * nCols + j];
        }
    }
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::set(uint r, uint c, U x){
    // allow to use -1 to indicate end of row/column
    if(r == UINT_MAX) r = this->_r-1;
    if(c == UINT_MAX) c = this->_c-1;

    this->at(r,c) = x;
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::set(uu_pair rs, uint c, Matrix<U> m){
    // allow to use -1 to indicate end of row/column
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    if(c == UINT_MAX) c = this->_c-1;

    if(rs.second >= this->_r) throw std::out_of_range("Row index out of range");
    if(c >= this->_c) throw std::out_of_range("Col index out of range");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second");

    if(m.shape() != uu_pair{rs.second - rs.first + 1, 1}) throw std::invalid_argument("Given matrix's shape does not match");

    for(uint i=0; i<m._r; ++i){
        this->at(i + rs.first, c) = m(i,0);
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::set(uint r, uu_pair cs, Matrix<T> m){
    // allow to use -1 to indicate end of row/column
    if(r == UINT_MAX) r = this->_r-1;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;

    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw std::out_of_range("Col index out of range");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(m.shape() != uu_pair{1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");

    for(uint j=0; j<m._c; ++j){
        std::copy(m._v, m._v + m.size(), this->_v + _c*r + cs.first);
    }
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::set(uint r, uu_pair cs, Matrix<U> m){
    // allow to use -1 to indicate end of row/column
    if(r == UINT_MAX) r = this->_r-1;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;

    if(r >= this->_r) throw std::out_of_range("Row index out of range");
    if(cs.second >= this->_c) throw std::out_of_range("Col index out of range");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(m.shape() != uu_pair{1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");

    for(uint j=0; j<m.c(); ++j){
        this->at(r, j + cs.first) = m(r,j);
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::set(uu_pair rs, uu_pair cs, Matrix<T> m){
    // allow to use -1 to indicate end of row/column
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;

    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(m.shape() != uu_pair{rs.second - rs.first + 1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");;

    for (uint i = 0; i < m.r(); ++i){
        std::copy(m._v + i * m.c(), m._v + (i+1) * m.c(), this->_v + (rs.first + i) * _c + cs.first);
    }
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::set(uu_pair rs, uu_pair cs, Matrix<U> m){
    // allow to use -1 to indicate end of row/column
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;

    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(m.shape() != uu_pair{rs.second - rs.first + 1, cs.second - cs.first + 1}) throw std::invalid_argument("Given matrix's shape does not match");;

    for (uint i = 0; i < m.r(); ++i){
        for (uint j = 0; j < m.c(); ++j){
            this->at(i + rs.first, j + cs.first) = m(i,j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::reshape(uint r, uint c) const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    if(r*c != this->size()) throw std::invalid_argument("New matrix size must match the current one");

    Matrix ret(r,c);
    std::copy(_v, _v + this->size(), ret._v);
    return ret;
}

template<typename T>
Matrix<T>& Matrix<T>::reshape_self(uint r, uint c){
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    if(r*c != this->size()) throw std::invalid_argument("New matrix size must match the current one");

    this->_r = r;
    this->_c = c;
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::diag(std::vector<U> v){
    if(std::min(_r,_c) != v.size()) throw std::invalid_argument("Size of v must match the size of the matrix diagonal");

    for(uint i=0; i<v.size(); ++i) this->at(i) = v[i];
    return *this;
}

template<typename T>
template<typename U, typename>
Matrix<T>& Matrix<T>::diag(Matrix<U> m){
    if(!m.is_vec()) throw std::invalid_argument("Input matrix must be vector-shaped");
    if(std::min(_r,_c) != m.size()) throw std::invalid_argument("Size of m must match the size of the matrix diagonal");

    for(uint i=0; i<m.size(); ++i) this->at(i) = m(i);
    return *this;
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
Matrix<T>& Matrix<T>::swap_rows(uint r1, uint r2){
    // allow to use -1 to indicate end of row/column
    if(r1 == UINT_MAX) r1 = this->_r-1;
    if(r2 == UINT_MAX) r2 = this->_r-1;

    if(r1 >= _r || r2 >= _r) throw std::out_of_range("Given parameters exceed the matrix rows indeces");

    // swap not necessary
    if(r1 == r2) return *this;

    // extract r1
    Matrix tmp = this->at(r1, ALL);
    // substitute r2 into r1
    this->set(r1, ALL, this->at(r2, ALL));
    // substitute r1 into r2
    this->set(r2, ALL, tmp);
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::swap_cols(uint c1, uint c2){
    // allow to use -1 to indicate end of row/column
    if(c1 == UINT_MAX) c1 = this->_c-1;
    if(c2 == UINT_MAX) c2 = this->_c-1;

    if(c1 >= _c || c2 >= _c) throw std::out_of_range("Given parameters exceed the matrix columns indeces");

    // swap not necessary
    if(c1 == c2) return *this;

    // extract c1
    Matrix tmp = this->at(ALL, c1);
    // substitute c2 into c1
    this->set(ALL, c1, this->at(ALL, c2));
    // substitute c1 into c2
    this->set(ALL, c2, tmp);
    return *this;
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
RetType_t<U,V> dot(const Matrix<U> & v1, const Matrix<V> & v2){
    if(! v1.is_vec()) throw std::invalid_argument("Obj 1 is not a vector");
    if(! v2.is_vec()) throw std::invalid_argument("Obj 2 is not a vector");
    if(v1.size() != v2.size()) throw std::invalid_argument("Vectors length don't match");

    RetType_t<U,V> ret = 0;
    for(uint i=0; i<v1.size(); ++i) ret += v1(i) * v2(i);
    return ret;
}

template<typename U, typename V>
Matrix<RetType_t<U,V>> cross(const Matrix<U> & v1, const Matrix<V> & v2){
    if(! v1.is_vec()) throw std::invalid_argument("Obj 1 is not a vector");
    if(! v2.is_vec()) throw std::invalid_argument("Obj 2 is not a vector");
    if(v1.size() != 3 || v1.size() != v2.size()) throw std::invalid_argument("Cross product defined only for 3d vectors");

    Matrix<RetType_t<U,V>> ret(3,1);
    ret(0) =   (v1(1) * v2(2)) - (v1(2) * v2(1));
    ret(1) = -((v1(0) * v2(2)) - (v1(2) * v2(0)));
    ret(2) =   (v1(0) * v2(1)) - (v1(1) * v2(0));
    return ret;
}

template<typename T>
double Matrix<T>::norm() const{
    if(!this->is_vec()) throw std::invalid_argument("Norm only appliable to row/column vectors");
    return std::sqrt(this->norm2());
}

template<typename T>
double Matrix<T>::norm2() const{
    if(!this->is_vec()) throw std::invalid_argument("Norm only appliable to row/column vectors");

    double ret = 0;
    for(uint i=0; i< this->size(); ++i){
        if constexpr(is_complex<T>::value) ret += _v[i].real() * _v[i].real() + _v[i].imag() * _v[i].imag();
        else ret += _v[i] * _v[i];
    }
    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::normalize() const{
    return Matrix<T>(*this) / this->norm();
}

template<typename T>
Matrix<T>& Matrix<T>::normalize_self(){
    return this->operator/=(this->norm());
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
    return abs(this->det()) < Matrix<T>::epsilon;
}

template<typename T>
bool Matrix<T>::is_upper_triang() const{
    for(uint i=1; i<_r; ++i){
        for(uint j=0; j<std::min<uint>(i, _c); ++j){
            if(abs(this->at(i,j)) > Matrix<T>::epsilon) return false;
        }
    }
    return true;
}
    
template<typename T>
bool Matrix<T>::is_lower_triang() const{
    for(uint j=1; j<_c; ++j){
        for(uint i=0; i<std::min<uint>(j, _r); ++i){
            if(abs(this->at(i,j)) > Matrix<T>::epsilon) return false;
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::is_upper_bidiagonal() const{
    if(!this->is_upper_triang()) return false;

    for(uint j=2; j<_c; ++j){
        for(uint i=0; i<std::min<uint>(j-1, _r); ++j){
            if(abs(this->at(i,j)) > Matrix<T>::epsilon) return false;
        }
    }

    return true;
}

template<typename T>
bool Matrix<T>::is_lower_bidiagonal() const{
    if(!this->is_lower_triang()) return false;

    for(uint i=2; i<_r; ++i){
        for(uint j=0; j<std::min<uint>(i-1, _c); ++j){
            if(abs(this->at(i,j)) > Matrix<T>::epsilon) return false;
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::is_upper_hessenberg() const{
    for(uint i=2; i<std::min(_r, _c); ++i){
        for(uint j=0; j<i-1; ++j){
            if(abs(this->at(i,j)) > Matrix<T>::epsilon) return false;
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::is_lower_hessenberg() const{
    for(uint j=2; j<std::min(_r, _c); ++j){
        for(uint i=0; i<j-1; ++i){
            if(abs(this->at(i,j)) > Matrix<T>::epsilon) return false;
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::is_orthogonal() const{
    return (*this) * this->t() == IdMat(_r);
}

#pragma endregion checks


#pragma region matrix_operations

template<typename T>
Matrix<T> Matrix<T>::t() const{
    Matrix<T> ret(this->_c, this->_r);
    for (uint i = 0; i<this->_r; ++i) {
        for (uint j = 0; j<this->_c; ++j) {
            if constexpr (is_complex<T>::value) ret._v[j * ret._c + i] = std::conj(this->_v[i * this->_c + j]);
            else ret._v[j * ret._c + i] = this->_v[i * this->_c + j];
        }
    }
    return ret;
}

template<typename T>
Matrix<T> Matrix<T>::no_conj_t() const{
    Matrix<T> ret(this->_c, this->_r);
    for (uint i = 0; i<this->_r; ++i) {
        for (uint j = 0; j<this->_c; ++j) {
            ret._v[j * ret._c + i] = this->_v[i * this->_c + j];
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
                ret(i,j) = this->at(row,col);
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

template<typename T>
Matrix<T> Matrix<T>::inv() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    return matrix_l_divide(*this, IdMat(_r));
}

template<typename T>
Matrix<T> Matrix<T>::pinv_left() const{
    if(_c > _r) throw std::invalid_argument("Matrix must be full column rank");
    try{
        Matrix<T> tmp = this->t() * *this;
        T det = tmp.det();
        if(abs(det) < Matrix<T>::epsilon) {
            std::cout << "Matrix is bad conditioned, det = " << det << "; this may result in bad numerical result" << std::endl;
        }
        return tmp.inv() * this->t();
    }
    catch(const std::runtime_error& rte){
        std::cerr << rte.what() << '\n';
        throw std::runtime_error("Matrix (m.t * m) not invertible");
    }
    catch(const std::exception& e){
        std::cerr << e.what() << '\n';
        throw std::runtime_error("Unknown error occured");
    }
}

template<typename T>
Matrix<T> Matrix<T>::pinv_right() const{
    if(_r > _c) throw std::invalid_argument("Matrix must be full row rank");
    try{
        Matrix<T> tmp = (*this * this->t());
        T det = tmp.det();
        if(abs(det) < Matrix<T>::epsilon) {
            std::cout << "Matrix is bad conditioned, det = " << det << "; this may result in bad numerical result" << std::endl;
        }
        return this->t() * tmp.inv();
    }
    catch(const std::runtime_error& rte){
        std::cerr << rte.what() << '\n';
        throw std::runtime_error("Matrix (m.t * m) not invertible");
    }
    catch(const std::exception& e){
        std::cerr << e.what() << '\n';
        throw std::runtime_error("Unknown error occured");
    }
}

template<typename T>
Matrix<c_double> Matrix<T>::sqrt() const{
    Matrix<c_double> V, D;
    this->eigen_dec(D, V);
    for(uint i=0; i<D.size(); ++i) D(i) = std::sqrt(D(i));    
    return V * MA::diag(D) * V.inv();
}

template<typename T>
Matrix<T> Matrix<T>::pow(int exp) const{
    if(_c != _r) throw std::invalid_argument("Matrix must be square");

    Matrix<T> z;
    Matrix<T> result = IdMat(_c);
    int bit;

    if(exp == 0) return IdMat(_c);
    else if(exp < 0){ 
        z = this->inv();
        exp = std::abs(exp);
    }
    else{
        z = *this;
    }

    while(exp > 0){
        bit = exp%2;
        exp = exp/2;

        if(bit) result *= z;

        z *= z;
    }

    return result;
}

#pragma endregion matrix_operations


#pragma region decomposition_methods

template<typename T>
Matrix<T> Matrix<T>::reflector(uu_pair rs, uint c) const{
    Matrix<T> u = this->at(rs, c);
    if(!u.is_vec()) throw std::invalid_argument("Given ranges do not extract a vector");

    using namespace std::complex_literals;
    T alpha;

    // if constexpr (is_complex<T>::value) alpha = -exp(arg(u(0))*1i) * u.norm();
    if constexpr (is_complex<T>::value) alpha = -u(0) / abs(u(0)) * u.norm();
    else alpha = -copysign(u.norm(), u(0));

    u(0) -= alpha;
    return u.normalize();
}

template<typename T>
Matrix<T> Matrix<T>::reflector(uint r, uu_pair cs) const{
    Matrix<T> u = this->at(r, cs).t();
    if(!u.is_vec()) throw std::invalid_argument("Given ranges do not extract a vector");

    using namespace std::complex_literals;
    T alpha;

    // if constexpr (is_complex<T>::value) alpha = -exp(arg(u(0))*1i) * u.norm();
    if constexpr (is_complex<T>::value) alpha = -u(0) / abs(u(0)) * u.norm();
    else alpha = -copysign(u.norm(), u(0));

    u(0) -= alpha;
    return u.normalize();
}

template<typename T>
void Matrix<T>::apply_reflector_right(
    const Matrix<T> & v, 
    uint start_col, 
    uu_pair rs
){
    // setup indexes
    uint i = start_col;
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    // range checks
    if(i + v.size() - 1 >= _c) throw std::invalid_argument("Reflector too big to be applied using given start_index");
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second");
    // init variables
    T tmp;
    Matrix<T> v_t = 2 * v.t();
    // apply reflector
    for(uint j=rs.first; j<=rs.second; ++j){
        tmp = 0.0;
        for(uint k=0; k<v.size(); ++k) tmp += v(k) * this->at(j,i+k);
        for(uint k=0; k<v.size(); ++k) this->at(j,i+k) -= v_t(k) * tmp;
    }
}

template<typename T>
void Matrix<T>::apply_reflector_left(
    const Matrix<T> & v, 
    uint start_row, 
    uu_pair cs
){
    // setup indexes
    uint i = start_row;
    if(cs.second == UINT_MAX) cs.second = this->_c-1;
    // range checks
    if(i + v.size() - 1 >= _r) throw std::invalid_argument("Reflector too big to be applied using given start_index");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");
    // init variables
    T tmp;
    Matrix<T> v_t = 2 * v.t();
    // apply reflector
    for(uint j=cs.first; j<=cs.second; ++j){
        tmp = 0.0;
        for(uint k=0; k<v.size(); ++k) tmp += v_t(k) * this->at(i+k,j);
        for(uint k=0; k<v.size(); ++k) this->at(i+k,j) -= v(k) * tmp;
    }
}

template<typename T>
Matrix<T> Matrix<T>::givens_rot(uint i, uint j) const{
    if(!this->is_vec()) throw std::invalid_argument("Given matrix must be vector shaped");
    Matrix<T> ret(2,1);
    T f = this->at(i);
    T g = this->at(j);

    if(abs(g) < Matrix<T>::epsilon && abs(f) < Matrix<T>::epsilon){
        ret(0) = 1.0;
        ret(1) = 0.0;
    }
    else{
        if constexpr (is_complex<T>::value){
            // T f2 = f.real()*f.real() + f.imag()*f.imag();
            // T d1 = 1.0/sqrt(f2 * (f2 + g.real()*g.real() + g.imag()*g.imag()));
            // ret(0) = f2 * d1;
            // ret(1) = f * d1 * std::conj(g);
            T d = 1/std::sqrt(f.real()*f.real() + f.imag()*f.imag() + g.real()*g.real() + g.imag()*g.imag());
            ret(0) = f * d;
            ret(1) = std::conj(g) * d;
        }
        else{
            T f2 = f*f;
            T d1 = 1.0/std::sqrt(f2 * (f2 + g*g));
            ret(0) = f2 * d1;
            ret(1) = f * d1 * g;
        }
    }
    return ret;
}

template<typename T>
void Matrix<T>::apply_givens_rot_right(const Matrix<T> & rot, uint i, uint j, uint r){
    if(!rot.is_vec() || rot.size() != 2) throw std::invalid_argument("Given rotation does not respect shape requirements");
    T f = this->at(r,i);
    T g = this->at(r,j);

    if constexpr(is_complex<T>::value) {
        this->at(r,i) = std::conj(rot(0))*f + rot(1)*g;
        this->at(r,j) = rot(0)*g - std::conj(rot(1))*f;
    }
    else {
        this->at(r,i) = rot(0)*f + rot(1)*g;
        this->at(r,j) = rot(0)*g - rot(1)*f;
    }
}

template<typename T>
void Matrix<T>::apply_givens_rot_right(const Matrix<T> & rot, uint i, uint j, uu_pair rs){
    if(rs.second == UINT_MAX) rs.second = this->_r-1;
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second");

    for(uint k=rs.first; k<=rs.second; ++k){
        this->apply_givens_rot_right(rot, i, j, k);
    }
}

template<typename T>
void Matrix<T>::apply_givens_rot_left(const Matrix<T> & rot, uint i, uint j, uint c){
    if(!rot.is_vec() || rot.size() != 2) throw std::invalid_argument("Given rotation does not respect shape requirements");
    T f = this->at(i,c);
    T g = this->at(j,c);
    
    if constexpr(is_complex<T>::value) {
        this->at(i,c) = std::conj(rot(0))*f + rot(1)*g;
        this->at(j,c) = rot(0)*g - std::conj(rot(1))*f;
    }
    else {
        this->at(i,c) = rot(0)*f + rot(1)*g;
        this->at(j,c) = rot(0)*g - rot(1)*f;
    }
}

template<typename T>
void Matrix<T>::apply_givens_rot_left(const Matrix<T> & rot, uint i, uint j, uu_pair cs){
    if(cs.second == UINT_MAX) cs.second = this->_c-1;
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    for(uint k=cs.first; k<=cs.second; ++k){
        this->apply_givens_rot_left(rot, i, j ,k);
    }
}

template<typename T>
void Matrix<T>::qr_dec(Matrix<T> & Q, Matrix<T> & R) const{
    // init matrices
    uint n = std::min<double>(_r, _c);
    Q = IdMat(_r);
    R = *this;
    
    for(uint i=0; i<n; ++i){
        //compute reflector
        Matrix<T> v = R.reflector({i, _r-1}, i); // U
        // apply reflector
        R.apply_reflector_left(v, i, {i,-1}); // UR
        Q.apply_reflector_right(v, i, ALL); // QU
    }
}

template<typename T>
void Matrix<T>::qrp_dec(Matrix<T> & Q, Matrix<T> & R, Matrix<double> & P) const{

    uint n = std::min<double>(_r, _c);

    // init matrices
    Q = IdMat(_r);
    R = *this;
    P = IdMat(_c);
    
    for(uint i=0; i<n; ++i){ // main loop

        // find column with largest norm
        uint index = 0;
        double max_norm = -1;
        for(uint k=i; k<n; ++k){
            double norm = R({i, _r-1}, k).norm();
            if(norm > max_norm){
                max_norm = norm;
                index = k;
            }
        }

        // swap columns
        R.swap_cols(i, index);
        P.swap_cols(i, index);

        // compute reflector
        Matrix<T> v = R.reflector({i, _r-1}, i); // U
        // apply reflector
        R.apply_reflector_left(v, i, {i,-1}); // UR
        Q.apply_reflector_right(v, i); // QU
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

    for (uint i=1; i<_r-1; ++i){
        // compute reflector
        Matrix<T> v = H.reflector({i, _r-1}, i-1); // U
        // apply reflector
        H.apply_reflector_left(v, i, {i-1,-1}); // UH
        H.apply_reflector_right(v, i); // HU
        Q.apply_reflector_right(v, i); // QU
    }
}

template<typename T>
void Matrix<T>::bidiagonal_form(Matrix<T> & U, Matrix<double> & B, Matrix<T> & Vt) const{
    if(_r < _c) throw std::invalid_argument("Number of rows must be greater or equal to number of columns");
    Matrix<T> B_tmp = *this;
    U = IdMat(_r,_r);
    Vt = IdMat(_c,_c);
    B = Matrix<double>(_r, _c);
    Matrix<T> v;

    for(uint i=0; i<_c; ++i){
        v = B_tmp.reflector({i, -1}, i); // Q
        B_tmp.apply_reflector_left(v, i, {i,-1}); // QB
        U.apply_reflector_right(v, i); // UQ
        if(i < _c-2){
            // need to transpose the vector becaue I am taking a row vector 
            // actually, it is necessary only for complex matrices
            v = B_tmp.reflector(i, {i+1, -1}); // Q
            B_tmp.apply_reflector_right(v, i+1, {i,-1}); // BQ
            Vt.apply_reflector_left(v, i+1); // QVt
        }
    }

    // if the matrix is complex, make the bidiagonal matrix real
    if constexpr (is_complex<T>::value){    
        for(uint i=0; i<_c; ++i){ // row index
            for(uint j=i; j<std::min(i+2, _c); ++j){ // column index
                double magn_inv = 1.0 / abs(B_tmp(i,j));
                T x(B_tmp(i,j).real() * std::copysign(magn_inv, B_tmp(i,j).imag()), - abs(B_tmp(i,j).imag()) * magn_inv);
                T x_c = std::conj(x);
                B(i,j) = (B_tmp(i,j) * x).real();
                if(i==j){ // modify row
                    // multiply row of B_tmp
                    if(j+1 < _c) B_tmp(i,j+1) *= x;
                    // multiply column of U
                    for(uint k=0; k<U.r(); ++k) U(k,j) *= x_c;
                }
                else{ // modify column
                    // multiply column of B_tmp
                    B_tmp(i+1,j) *= x;
                    // multiply row of Vt
                    for(uint k=0; k<Vt.c(); ++k) Vt(j,k) *= x_c;
                }
            }
        }
        // std::cout << B << std::endl;
    }
    else{
        B = B_tmp;
    }
}

template<>
double Matrix<double>::svd_shift() const{
    uint n = std::min(_r-1,_c-1);
    using namespace std;

    double a, b, c, d;
    // matrix elements: | a  b |
    //                  | c  d |
    a = this->at(n-1)*this->at(n-1) + this->at(n-2,n-1)*this->at(n-2,n-1);
    b = this->at(n-1) * this->at(n-1,n);
    c = b;
    d = this->at(n)*this->at(n) + this->at(n-1,n)*this->at(n-1,n);
    // polynom coefficients
    c = a * d - c * b;
    b = (- a - d)/2.0; // b/2
    // matrix is symmetric -> only real eigenvalues 
    a = std::sqrt(b*b - c); // alpha
    c = -b+a;
    d = -b-a;
    // return closest to A(n,n)
    if(abs(c - this->at(n)) <= abs(d - this->at(n))) return c;
    else return d;
}

template<typename T>
void Matrix<T>::svd_reinsch_step(
    Matrix<T> & P,
    Matrix<double> & E,
    Matrix<T> & Gt,
    double shift,
    uu_pair range
){
    if(E.r() < E.c()) throw std::invalid_argument("Matrix must represent an overdetermined problem");
    if(range.second == UINT_MAX) range.second = E.c()-1;
    if(range.second >= E.c()) throw std::out_of_range("Given bounds exceed matrix dimensions");
    if(range.first > range.second) throw std::invalid_argument("Range first element must be <= of second"); 
    uint start = range.first;
    uint end = range.second;

    using namespace std;

    // apply shift
    // compute first reflector
    Matrix<double> v = Matrix<double>(2, 1, {E(start)*E(start) - shift, E(start)*E(start,start+1)}).givens_rot(0,1); // G1
    E.apply_givens_rot_right(v, start, start+1, {start, std::min(start+2, end)}); // EG1
    Gt.apply_givens_rot_left(v, start, start+1); // G1Gt

    // cout << "after shift: " << endl << E << endl;
    // chase the bulge
    for(uint i=start; i<end; ++i){
        v = E({i, i+1}, i).givens_rot(0,1); // Pi
        E.apply_givens_rot_left(v, i, i+1, {i, std::min(i+2, end)}); // PiE
        P.apply_givens_rot_right(v, i, i+1); // PPi
        // cout << "after lower bulge: " << endl << E << endl;

        if(i < end-1){
            v = E(i, {i+1, i+2}).givens_rot(0,1); // Gi
            E.apply_givens_rot_right(v, i+1, i+2, {i, std::min(i+2, end)}); // EGi
            Gt.apply_givens_rot_left(v, i+1, i+2); // GiGt
            // cout << "after upper bulge: " << endl << E << endl;
        }
    }

}

template<typename T>
void Matrix<T>::svd_steps_iteration(
    Matrix<T> & P,
    Matrix<double> & E,
    Matrix<T> & Gt,
    uu_pair range,
    uint max_iterations,
    double tolerance
){
    if(E.r() < E.c()) throw std::invalid_argument("Matrix must represent an overdetermined problem");
    if(range.second == UINT_MAX) range.second = E.c()-1;
    if(range.second >= E.c()) throw std::out_of_range("Given bounds exceed matrix dimensions");
    if(range.first > range.second) throw std::invalid_argument("Range first element must be <= of second"); 
    // uint si = start_index;
    
    uint zero_shift_iterations = 0;
    double shift;
    bool converged = false;

    using namespace std;

    // cout << "starting index: " << si << ", dim: " << dim << endl
    //      << "working matrix:" << E({si,si+dim-1},{si,si+dim-1}) << endl;

    // I can stop: we deflated a single value, meaning it converged
    if(range.first == range.second){ 
        // cout << "single value: iteration terminated" << endl;
        return;
    }

    for(uint steps=0; steps<max_iterations && !converged; ++steps){

        // deflation
        for(uint i=range.second; i>=range.first+1; --i){ // i indicates column position
            if(abs(E(i-1,i)) < tolerance){
                // cout << P * E * Gt << endl << "converged in i=" << i << endl << E << endl;
                converged = true;
                // cout << steps << endl;
                // upper part
                Matrix<T>::svd_steps_iteration(P, E, Gt, {range.first, i-1}, max_iterations, tolerance);
                // lower part
                Matrix<T>::svd_steps_iteration(P, E, Gt, {i, range.second}, max_iterations, tolerance);
                // exit
                return;
            }
        }
        // cancellation
        // for(uint i=range.second; i>=range.first+1; --i){
        //     // ATT: without +1 it cycle for inf (once it's 0, it overflows)
        // }

        // apply step
        if(steps < zero_shift_iterations) shift = 0;
        else shift = E(range, range).svd_shift();
        // cout << "matrix before step: " <<  endl << E << endl;
        Matrix<T>::svd_reinsch_step(P, E, Gt, shift, range);
    }
}

template<typename T>
void Matrix<T>::svd(Matrix<T> & U, Matrix<double> & E, Matrix<T> & Vt, uint max_iterations, double tolerance) const{
    using namespace std;

    // transpose this if r<c
    Matrix<T> A;
    if(_r<_c) A = this->t();
    else A = *this;

    // reduce to bidiagonal form
    A.bidiagonal_form(U,E,Vt);
    // cout << "starting matrices: " << E 
    //      << endl << U << endl << Vt << endl;
    
    // init matrices
    Matrix<T> P = IdMat(E.r());
    Matrix<T> Gt = IdMat(E.c());

    Matrix<T>::svd_steps_iteration(P,E,Gt, ALL, max_iterations, tolerance);
    // cout << P << endl << E << endl << Gt << endl;
    // make all singular values have positive values
    for(uint i=0; i<E.c(); ++i){
        if(E(i) < 0){
            E(i) = -E(i);
            for(uint j=0; j<Gt.c(); ++j){
                Gt(i,j) = -Gt(i,j);
            }
        }
    }

    // update solution incorporating bidiagonal decomposition
    U = U*P;
    Vt = Gt*Vt;

    // reorder singular values
    // create vector of singular values + indeces
    std::vector<std::pair<double,uint>> vec_backward_sort;
    std::vector<uint> vec_forward_sort;
    for(uint i=0; i<E.c(); ++i){
        vec_backward_sort.push_back(std::pair<double, uint>(E(i), i));
        vec_forward_sort.push_back(0);
    }
    auto sort_logic = [](std::pair<double, uint> i, std::pair<double, uint> j) {
        return i.first > j.first;
    };
    // reorder the vector to keep track of the indeces
    std::sort(vec_backward_sort.begin(), vec_backward_sort.end(), sort_logic);
    // define the destination of each element from non-sorted to sorted
    for(uint i=0; i<vec_backward_sort.size(); ++i){
        vec_forward_sort[vec_backward_sort[i].second] = i;
    }
    // reorder matrices
    uint i=0;
    while(i<vec_forward_sort.size()){
        if(vec_forward_sort[i] == i) ++i;
        else{
            U.swap_cols(i, vec_forward_sort[i]);
            E.swap_rows(i, vec_forward_sort[i]);
            E.swap_cols(i, vec_forward_sort[i]);
            Vt.swap_rows(i, vec_forward_sort[i]);
            std::swap(vec_forward_sort[i], vec_forward_sort[vec_forward_sort[i]]);
        }
    }

    // transpose the result
    if(_r<_c){
        E = E.t();
        Matrix<T> tmp = U;
        U = Vt.t();
        Vt = tmp.t();
    }
}

#pragma endregion decomposition_methods


#pragma region ls_solution

template<typename T, typename V, typename W>
Matrix<W> backward_sub(Matrix<T> const & U, Matrix<V> const & B){
    if(U.c() != U.r()) throw std::invalid_argument("Coefficient matrix U must be square");
    if(B.r() != U.c()) throw std::invalid_argument("Rows of B must be equal to the columns of U");
    // check all U diagonal elements are different from zero
    for(uint i=0; i<U.c(); ++i) if(U(i,i) == 0.0){
        throw std::runtime_error("System of equation is underdetermined");
    }

    Matrix<W> res(U.r(), B.c());

    W tmp;
    // printf("start i loop\n");
    for(uint i=0; i<B.c(); ++i){ //col of res == col of B
        // printf("start j loop with i:%d\n", i);
        for(int j=U.r()-1; j>=0; --j){ //row of res == row of U == row of B
            tmp = 0.0;
            // printf("start k loop with i:%d, j:%d\n", i, j);
            for(int k=U.r()-1; k>j; --k){ //col of U = row of res
                tmp += U(j,k) * res(k,i);
                // printf("i:%d j:%d k:%d\n",i,j,k);
            }
            
            res(j, i) = (B(j,i) - tmp) / U(j,j);
            // printf("%f\t%f\n",tmp,res[j*B.c() + i]);
        }
        // printf("\n");
    }
    
    return res;
}

template<typename T, typename U, typename V>
Matrix<V> forward_sub(Matrix<T> const & L, Matrix<U> const & B){
    if(L.c() != L.r()) throw std::invalid_argument("Coefficient matrix L must be square");
    if(B.r() != L.c()) throw std::invalid_argument("Rows of B must be equal to the columns of L");
    // check all L diagonal elements are different from zero
    for(uint i=0; i<L.c(); ++i) if(L(i,i) == 0.0){
        throw std::runtime_error("System of equation is underdetermined");
    }

    Matrix<V> res(L.r(), B.c());
    
    V tmp;
    for(uint j=0; j<L.r(); ++j){ //row of res == row of L == row of B
        // printf("j: %d\n", j);
        for(uint i=0; i<B.c(); ++i){ //col of res == col of B
            // printf("i: %d\n", i);
            tmp = 0.0;
            for(uint k=0; k<j; ++k){ //col of L = row of res
                // printf("k: %d\n", k);
                tmp += L(j,k) * res(k,i);
            }
            // printf("%f\n", tmp);
            res(j,i) = (B(j,i) - tmp) / L(j,j);
        }
    }
    
    return res;
}

template<typename T, typename U, typename V>
Matrix<V> matrix_l_divide(Matrix<T> const & A, Matrix<U> const & B){
    if(A.r() != B.r()){
        throw std::invalid_argument("Rows of B (currently " + std::to_string(B.r()) + ")" +
        "must be equal the rows of A (currently " + std::to_string(A.r()) + ")");
    }

    if(A.c() == A.r()){ //square A
        // printf("matrix A is square, using LU decomposition\n");
        Matrix<T> L, u;
        Matrix<double> P;
        Matrix<V> res;

        A.lup_dec(L, u, P);
        res = forward_sub(L, P*B);
        return backward_sub(u, res);
    }
    else{
        // printf("matrix A is not square, using QR decomposition\n");
        if(A.r() < A.c()) throw std::invalid_argument("System is underdetermined");
        Matrix<T> Q, R;
        Matrix<double> P;
        Matrix<V> res;

        A.qrp_dec(Q,R,P);
        res = P * backward_sub(R({0, R.c()-1}, ALL), (Q.t() * B)({0, R.c()-1}, ALL));
        return res;
    }
}

template<typename T, typename U, typename V>
Matrix<V> solve_ls(Matrix<T> const & A, Matrix<U> const & B){
    return matrix_l_divide(A,B);
}

template<typename T, typename U, typename V>
Matrix<V> matrix_r_divide(Matrix<T> const & B, Matrix<U> const & A){
    if(A.c() != B.c()){
        throw std::invalid_argument("Columns of B (currently " + std::to_string(B.c()) + ") " +
        "must be equal the Columns of A (currently " + std::to_string(A.c()) + ")");
    }
    return matrix_l_divide(A.no_conj_t(), B.no_conj_t()).no_conj_t();
}

#pragma endregion ls_solution


#pragma region eigen

template<typename T>
Matrix<T> Matrix<T>::implicit_double_QR_step() const{
    // if(!this->is_upper_hessenberg()) throw std::invalid_argument("Input matrix must be upper hessenberg");

    // for simpler notation
    uint n =  this->_r;
    uint start;
    Matrix<T> A = *this;

    Matrix<T> y, v;

    for(uint i=0; i<n-1; ++i){
        if(i == 0){
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
            start = 0;
        }
        else {
            // extract column to transform in [x 0 0]

            y = A({i, std::min(i+2, n-1)}, i-1);
            start = i-1;
        }
        // compute reflector
        v = y.reflector(); // Q
        // apply reflector
        A.apply_reflector_left(v, i, {start, -1});
        A.apply_reflector_right(v, i);
    }

    return A;
}

template<typename T>
Matrix<c_double> Matrix<T>::eigenvalues(uint max_iterations, double tolerance) const{
    if(_c != _r) throw std::invalid_argument("Matrix must be square");

    using namespace std::complex_literals;

    // if(_r == 0) throw std::invalid_argument("Matrix is empty");
    if(_r <= 1) return *this;
    else if(_r == 2){
        T b = -this->at(0,0) - this->at(1,1);
        T c = this->at(0,0) * this->at(1,1) - this->at(0,1) * this->at(1,0);
        T alpha = b*b - 4.0*c;
        
        if constexpr (is_complex<T>::value){
            T tmp = std::sqrt(alpha);
            return Matrix<c_double>(2,1, std::vector<c_double>{
                (- b + tmp) / 2.0,
                (- b - tmp) / 2.0
            });
        }
        else{
            T tmp = std::sqrt(abs(alpha));
            if(alpha >= 0.0){
                return Matrix<c_double>(2,1, std::vector<c_double>{
                    (- b + tmp) / 2.0,
                    (- b - tmp) / 2.0
                });
            }
            else{
                return Matrix<c_double>(2,1, std::vector<c_double>{
                    (- b + tmp*1i) / 2.0,
                    (- b - tmp*1i) / 2.0
                });
            }
        }
    }
    else{
        Matrix<T> Q, H;
        // obtain an hessenberg matrix
        if(!this->is_upper_hessenberg()) {
            this->hessenberg_dec(Q,H);
        }
        else{
            H = *this;
            Q = IdMat(_r);
        }

        // start running the implicit double QR steps
        for(uint k=0; k<max_iterations; ++k){
            H = H.implicit_double_QR_step();
            // check for convergence
            for(int i=_c-2; i>=0; --i){
                if(abs(H(i+1,i)) < tolerance){
                    // deflation
                    Matrix<c_double> m1 = H(uu_pair{0,i},uu_pair{0,i});
                    Matrix<c_double> m2 = H(uu_pair{i+1,_r-1},uu_pair{i+1,_r-1});
                    return m1.eigenvalues(max_iterations, tolerance) | m2.eigenvalues(max_iterations, tolerance);
                }
            }
        }

        std::cout << "eigenvalue extraction failed to converge" << std::endl
                  << H << std::endl;
        return H.diag();
    }
}

template<typename T>
Matrix<c_double> Matrix<T>::eigenvector(c_double eigenv, uint max_iterations, double tolerance) const{
    if(_c != _r) throw std::invalid_argument("Matrix must be square");

    /*
    solving the linear system associated the equation A*v=l*v
    prooves to be as precise as shifted inverse iteration.
    However the reasoning is that shifted inverse iteration converges
    to the eigenvector even if the eigenvalues is subject to small errors
    */
    // // extract eigenvector solving the linear system
    // Matrix<c_double> M = *this - eigenv * IdMat(_r);
    // Matrix<c_double> eigen_vec = Matrix(1,1,{1}) | solve_ls(M(ALL, {1,_r-1}), - M(ALL, 0));
    // // std::cout << (*this * eigen_vec - eigenv * eigen_vec).norm() << std::endl;
    // return eigen_vec;
    
    // compute LU decomposition of shifted matrix
    Matrix<c_double> L, U, tmp;
    Matrix<double> P;
    (*this - IdMat(_r) * eigenv).lup_dec(L,U,P);
    uint k;
    double residual, prev_residual = 1;

    // start inverse iteration
    Matrix<c_double> q = Ones(_r,1).normalize();
    // Matrix<c_double> q = RandMat<c_double>(_r,1);
    for(k=0; k<max_iterations; ++k){
        tmp = forward_sub(L, P*q);
        q = backward_sub(U, tmp).normalize();
        residual = (*this * q - eigenv * q).norm();
        // check for convergence
        if(residual < tolerance) break;
        // the algorithm should converge very fast
        if(prev_residual / residual < 10) break;

        prev_residual = residual;
    }
    
    // std::cout << k << "; " << residual << std::endl;
    return q;
}

template<typename T>
void Matrix<T>::eigen_dec(Matrix<c_double> & D, Matrix<c_double> & V, uint max_iterations, double tolerance) const{
    if(_c != _r) throw std::invalid_argument("Matrix must be square");

    // extract eigenvalues
    D = this->eigenvalues(max_iterations, tolerance);
    // extract eigenvectors using shifted inverse iteration
    V = Matrix<c_double>(_r, _c);
    for(uint i=0; i<D.size(); ++i){
        V.set(ALL, i, this->eigenvector(D(i)));
    }

}

#pragma endregion eigen


#pragma region special_constructors

Matrix<double> IdMat(uint r, uint c){
    Matrix ret = Matrix(r,c);
    for(uint i=0; i<std::min(r,c); ++i) ret(i) = 1;
    return ret;
}

Matrix<double> IdMat(uint dim){
    return IdMat(dim, dim);
}

Matrix<double> Ones(uint r, uint c){
    return Matrix(r,c)+1;
}

template<typename T>
Matrix<T> RandMat(uint r, uint c){
    Matrix<T> ret(r,c);
    for(uint i=0; i<r; ++i){
        for(uint j=0; j<c; ++j){
            ret(i,j) = Matrix<T>::rand();
        }
    }

    return ret;
}

template<typename T>
Matrix<T> diag(const std::vector<T> & v){
    Matrix<T> ret(v.size(), v.size());
    
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
