#ifndef MA_MATRICES_HPP
#define MA_MATRICES_HPP

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cmath>
#include <climits> 
#include <complex>
#include <stdexcept>
#include <exception>

#include <cstddef>
#include <concepts>
#include <type_traits>

#include <iostream>
#include <sstream>
#include <string>

#include <utility>
#include <vector>
#include <random>


namespace MA
{

typedef unsigned int uint;
typedef std::pair<uint,uint> uu_pair;
typedef std::complex<double> c_double;
// #define ALL uu_pair{0, UINT_MAX}
#define ALL uu_pair{}

#pragma region helper_structs

// Define a type trait to check if a type is std::complex
template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

// Helper type to conditionally determine the type of V
template <typename T, typename U>
struct RetType {
    using type = std::conditional_t<is_complex<T>::value || is_complex<U>::value, c_double, double>;
};
template <typename T, typename U>
using RetType_t = typename RetType<T, U>::type;

// Helper type to conditionally enable a function (if there is no conversion from complex to double)
template <typename T, typename U, typename = void>
struct enable_if_not_comp2d : 
    std::enable_if<
        !is_complex<U>::value || 
        (is_complex<T>::value && is_complex<U>::value),
    int> 
{};
template <typename T, typename U>
using enable_if_not_comp2d_t = typename enable_if_not_comp2d<T, U>::type;

#pragma endregion helper_structs

template<typename T = double>
class Matrix {
protected:

    // number of decimals used in double comparison
    inline static double double_precision = 16.0;
    // epsilon used in double comparison
    inline static double epsilon = 1e-16;
    // random number generator objects
    inline static std::uniform_real_distribution<double> unif{0.0, 1.0};
    inline static std::default_random_engine re;

    uint _r; // rows
    uint _c; // columns
    T * _v; // vector storing the elements per rows

public:

#pragma region static_methods
    /**
     * @brief set the double precision used in comparison
     * @param dp number of digits considered during comparison
     */
    inline static void set_double_precision(double dp = 16.0){ 
        if(dp > 16) std::cout << "WARNING: double precision very small, this may result in bad behavior" << std::endl;
        Matrix<T>::double_precision = dp; 
        Matrix<T>::epsilon = pow(10.0,-dp);
    }

    /**
     * @brief get the double precision used in comparison
     * @return margin of error considered during comparison
     */
    inline static double get_epsilon(){ 
        return Matrix<T>::epsilon; 
    }

    /**
     * @brief generate a random number using the internally set seed
     * @return random number in range 0..1
     */
    static T rand();
#pragma endregion static_methods
 
#pragma region constructor_destructor
    /**
     * @brief empty constructor, rows and columns are set to 0 
     *      and v is set to nullptr 
    */
    Matrix();

    /**
     * @brief construct a r x c matrix with all elements 
     *      initialized to zero.
     * @param r number of rows
     * @param c number of columns
    */
    Matrix(const uint & r, const uint & c);

    /**
     * @brief construct a r x c matrix using v to initialize
     *      the elements (per rows).
     * @param r number of rows
     * @param c number of columns
     * @param v vector used to init the matrix
     */
    Matrix(const uint & r, const uint & c, std::vector<T> v);

    /**
     * @brief construct a r x c matrix using v to initialize
     *      the elements (per rows).
     * @param r number of rows
     * @param c number of columns
     * @param v vector used to init the matrix
     */
    template<typename U,
        typename std::enable_if<
            !is_complex<U>::value ||
            (is_complex<T>::value && is_complex<U>::value), int
        >::type = 0
    >
    Matrix(const uint & r, const uint & c, std::vector<U> v);

    /**
     * @brief copy constructor
     * @param m matrix to copy
     */
    Matrix(const Matrix<T> & m);

    /**
     * @brief copy constructor
     * @param m matrix to copy
    */
    template<typename U,
        typename = enable_if_not_comp2d_t<T,U>
    >
    Matrix(const Matrix<U> & m);

    /**
     * @brief move constructor
    */
    Matrix(Matrix<T> && m) noexcept;

    /**
     * @brief destructor, makes sure to free the memory
     *      allocated for v
    */
    ~Matrix();
#pragma endregion constructor_destructor

#pragma region get
    /**
     * @brief get operator using ()
     * @param r row index
     * @param c columns index
     * @return T& element in position (r,c)
     */
    T& operator()(const uint & r, const uint & c);

    /**
     * @brief const get operator using ()
     * @param r row index
     * @param c columns index
     * @return const T& element in position (r,c)
     */
    const T& operator()(const uint & r, const uint & c) const;

    /**
     * @brief if the object is a vector, access the i-th element,
     * if the object is a matrix, access the i-th element on the diagonal.
     * The returned value is a modifiable reference.
     * @param i index of access
     * @return element of the object
    */
    T& operator()(const uint & i);

    /**
     * @brief if the object is a vector, access the i-th element,
     * if the object is a matrix, access the i-th element on the diagonal.
     * The returned value is a constant refrence.
     * @param i index of access
     * @return element of the object
    */
    const T& operator()(const uint & i) const;

    /**
     * @brief extract the elements on rows rs and column c.
     * @param rs rows index pair (both extremes included)
     * @param c column index
     * @return submatrix of rows rs and column c
    */
    Matrix<T> operator()(uu_pair rs, const uint & c) const;

    /**
     * @brief extract the elements on row r and columns cs.
     * @param r row index 
     * @param cs columns index pair (both extremes included)
     * @return submatrix of rows r and column cs
    */
    Matrix<T> operator()(const uint & r, uu_pair cs) const;

    /**
     * @brief extract the elements on rows rs and columns cs.
     * @param rs rows index pair (both extremes included)
     * @param cs columns index pair (both extremes included)
     * @return submatrix of rows rs and columns cs
    */
    Matrix<T> operator()(uu_pair rs, uu_pair cs) const;

    /**
     * @brief get number of rows
     * @return number of rows
     */
    uint r() const;

    /**
     * @brief get number of columns
     * @return number of columns
     */
    uint c() const;

    /**
     * @brief return matrix shape as (r,c) pair
     * @return pair
     */
    uu_pair shape() const;

    /**
     * @brief get size of the matrix i.e. number of stored elements
     * @return size
    */
    uint size() const;

    /**
     * @brief return underlaying array storing the elements
     * @return underlaying array
     */
    T const * v() const;

    /**
     * @brief extract diagonal of the matrix
     * @return vector (matrix with dim n,1) containing the diagonal elements
    */
    Matrix<T> diag() const;

    /**
     * @brief return matrix containing real part of the complex matrix
    */
    template <typename U = T, typename std::enable_if<is_complex<U>::value, int>::type = 0>
    Matrix<double> real() const;

    /**
     * @brief return matrix containing imaginary part of the complex matrix
    */
    template <typename U = T, typename std::enable_if<is_complex<U>::value, int>::type = 0>
    Matrix<double> imag() const;

#pragma endregion get

#pragma region set
    Matrix<T>& operator=(Matrix<T> m);

    /**
     * @brief uses the vector to fill in the elements of the matrix
     * @param v used vector
     * @throw out_of_range if v.size < this.size
     */
    Matrix<T>& set(std::vector<T> v);

    /**
     * @brief uses the vector to fill in the elements of the matrix
     * @param v used vector
     * @throw out_of_range if v.size < this.size
     */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& set(std::vector<U> v);

    /**
     * @brief uses the vector to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns vector
     * @param v used vector
     * @throw out_of_range if v.size < submatrix size
     */
    Matrix<T>& set(uu_pair rs, uu_pair cs, std::vector<T> v);

    /**
     * @brief uses the vector to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns vector
     * @param v used vector
     * @throw out_of_range if v.size < submatrix size
     */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& set(uu_pair rs, uu_pair cs, std::vector<U> v);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs
     * @param r row index
     * @param c column index
     * @param m used matrix
     * @throw out_of_range if r or c are out of range
     */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& set(uint r, uint c, U x);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces rs and c, extremes included
     * @param rs rows indeces
     * @param c column index
     * @param m used matrix
     * @throw out_of_range if rs or c are out of range
     * @throw out_of_range if m has not enough rows or columns to fill in the 
     * submatrix
     */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& set(uu_pair rs, uint c, Matrix<U> m);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces r and cs, extremes included
     * @param r row index
     * @param cs columns indeces
     * @param m used matrix
     * @throw out_of_range if r or cs are out of range
     * @throw out_of_range if m has not enough rows or columns to fill in the 
     * submatrix
     */
    Matrix<T>& set(uint r, uu_pair cs, Matrix<T> m);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces r and cs, extremes included
     * @param r row index
     * @param cs columns indeces
     * @param m used matrix
     * @throw out_of_range if r or cs are out of range
     * @throw out_of_range if m has not enough rows or columns to fill in the 
     * submatrix
     */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& set(uint r, uu_pair cs, Matrix<U> m);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns indeces
     * @param m used matrix
     * @throw out_of_range if rs or cs are out of range
     * @throw out_of_range if m has not enough rows or columns to fill in the 
     * submatrix
     */
    Matrix<T>& set(uu_pair rs, uu_pair cs, Matrix<T> m);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns indeces
     * @param m used matrix
     * @throw out_of_range if rs or cs are out of range
     * @throw out_of_range if m has not enough rows or columns to fill in the 
     * submatrix
     */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& set(uu_pair rs, uu_pair cs, Matrix<U> m);

    /**
     * @brief returns a matrix with the same underlaying
     * array, but different shape. Size must be the same
     * @param r rows of new matrix
     * @param c cols of new matrix
     * @return Matrix: matrix with same data and specified shape
     */
    Matrix<T> reshape(const uint & r, const uint & c) const;

    /**
     * @brief Modifies the shape of this matrix. Size must be the same
     * @param r rows of new matrix
     * @param c cols of new matrix
     * @return ref to this matrix
     */
    Matrix<T>& reshape_self(const uint & r, const uint & c);

    /**
     * @brief swap as friend function to allow ADL
     * @param m1 first matrix
     * @param m2 second matrix
     */
    template<typename U>
    friend void swap(Matrix<U> & m1, Matrix<U> & m2);

    /**
     * @brief swap the given rows of the matrix
     * @param r1 first row to swap
     * @param r2 second row to swap
     * @throw invalid argument if r1 or r2 exceeds the matric row indeces
     */
    Matrix<T>& swap_rows(const uint & r1, const uint & r2);

    /**
     * @brief swap the given columns of the matrix
     * @param c1 first column to swap
     * @param c2 second column to swap
     * @throw invalid argument if c1 or c2 exceeds the matric column indeces
     */
    Matrix<T>& swap_cols(const uint & c1, const uint & c2);
    
    /**
     * @brief set diagonal elements of this matrix using given vector
     * @param v vector used to set the diagonal
     * @return ref to this matrix
     * @throw if dim of v does not match the dim of the diagonal of the matrix
    */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& diag(std::vector<U> v);

    /**
     * @brief set diagonal elements of this matrix using given matrix
     * @param m matrix used to set the diagonal, must have a vector shape
     * @return ref to this matrix
     * @throw if dim of m does not match the dim of the diagonal of the matrix
    */
    template<typename U, typename = enable_if_not_comp2d_t<T,U>>
    Matrix<T>& diag(Matrix<U> m);
#pragma endregion set

#pragma region comparison_operators

    template <typename U, typename V>
    friend std::enable_if_t<!is_complex<U>::value && !is_complex<V>::value, bool>
    operator==(const Matrix<U> & m1, const Matrix<V> & m2);

    template <typename U, typename V>
    friend std::enable_if_t<is_complex<U>::value || is_complex<V>::value, bool>
    operator==(const Matrix<U> & m1, const Matrix<V> & m2);

    template<typename U, typename V>
    friend bool operator!=(const Matrix<U> & m1, const Matrix<V> & m2);

#pragma endregion comparison_operators

#pragma region sum_operator
    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator+=(const V & k);

    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator+=(const Matrix<V> & m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator+(const Matrix<U>& m, const V& k);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator+(const U& k, const Matrix<V>& m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator+(const Matrix<U>& m1, const Matrix<V>& m2);
#pragma endregion sum_operator

#pragma region subtract_operator
    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator-=(const V & k);

    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator-=(const Matrix<V> & m);

    template<typename U>
    friend Matrix<U> operator-(const Matrix<U>& m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator-(const Matrix<U>& m, const V& k);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator-(const U& k, const Matrix<V>& m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator-(const Matrix<U>& m1, const Matrix<V>& m2);
#pragma endregion subtract_operator

#pragma region multiply_operator
    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator*=(const V & k);

    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator*=(const Matrix<V> & m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator*(const Matrix<U>& m, const V& k);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator*(const U& k, const Matrix<V>& m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator*(const Matrix<U>& m1, const Matrix<V>& m2);
#pragma endregion multiply_operator

#pragma region divide_operator
    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator/=(const V & k);

    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator/=(const Matrix<V> & m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator/(const Matrix<U>& m, const V& k);
    
    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator/(const U& k, const Matrix<V>& m);

    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator/(const Matrix<U>& m1, const Matrix<V>& m2);
#pragma endregion divide_operator

#pragma region output_operator
    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& m);
#pragma endregion output_operator

#pragma region concatenate_operators
    /**
     * @brief concatenate matrices per columns
     * @param m matrix to concatenate
     */
    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator&=(const Matrix<V> & m);

    /**
     * @brief concatenate matrices per columns
     * @param m1 first matrix to concatenate
     * @param m2 second matrix to concatenate
     */
    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator&(const Matrix<U>& m1, const Matrix<V>& m2);

    /**
     * @brief concatenate matrices per rows
     * @param m matrix to concatenate
     */
    template <typename U = T, typename V, typename = enable_if_not_comp2d_t<U,V>>
    Matrix<T>& operator|=(const Matrix<V> & m);

    /**
     * @brief concatenate matrices per rows
     * @param m1 first matrix to concatenate
     * @param m2 second matrix to concatenate
     */
    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> operator|(const Matrix<U>& m1, const Matrix<V>& m2);
#pragma endregion concatenate_operators

#pragma region vector
    /**
     * @brief reshape matrix into column vec
     * @return Matrix: column vec
     */
    Matrix<T> to_c_vec() const;

    /**
     * @brief reshape matrix into row vec
     * @return Matrix: row vec
     */
    Matrix<T> to_r_vec() const;

    /**
     * @brief computes dot product this.v
     * @param v second vector
     * @return dot product 
     */
    template<typename U, typename V>
    friend RetType_t<U,V> dot(const Matrix<U> & v1, const Matrix<V> & v2);

    /**
     * @brief computes cross product this*v
     * implemented only for vectors of length 2 and 3
     * @param v second vector
     * @return Matrix: cross product (as column vec)
     */
    template<typename U, typename V>
    friend Matrix<RetType_t<U,V>> cross(const Matrix<U> & v1, const Matrix<V> & v2);

    /**
     * @brief return norm2: sqrt(sum(v(i)^2)) of the given vector
     * 
     * @return double: norm2 of the vector
     */
    double norm2() const;

    /**
     * @brief return the normalized vector
     * 
     * @return Matrix: normalized vector (matrix with one dim=1)
     */
    Matrix<T> normalize() const;

    /**
     * @brief normalized the vector, modifying current object
     * 
     */
    Matrix<T>& normalize_self();
#pragma endregion vector

#pragma region checks
    /**
     * @brief check if matrix is a vector
     * i.e. one of the dimension is 1
     * @return true if is at least one dim is == 1
     */
    bool is_vec() const;

    /**
     * @brief test if the given matrix is singular
     * matrix must be square
     * @return bool: true if it's singular
     */
    bool is_sing() const;

    /**
     * @brief return true if the matrix is upper triangular
     */
    bool is_upper_triang() const;
    
    /**
     * @brief return true if the matrix is lower triangular
     */
    bool is_lower_triang() const;

    /**
     * @brief return true if the matrix is an upper hessenberg matrix
     */
    bool is_upper_hessenberg() const;

    /**
     * @brief return true if the matrix is an lower hessenberg matrix
     */
    bool is_lower_hessenberg() const;
#pragma endregion checks

#pragma region matrix_operations
    /**
     * @brief compute transpose, if the matrix is complex
     *  it computes the conjugate transpose
     * @return Matrix: transpose
     */
    Matrix<T> t() const;

    /**
     * @brief compute simple transpose, even if the matrix is complex
     *  (does not conjugate the entries)
     * @return Matrix: transpose
     */
    Matrix<T> no_conj_t() const;

    /**
     * @brief compute submatrix matrix obtained by deleting
     * the p-th row and q-t column
     * 
     * @param p row index
     * @param q column index
     * @return Matrix: submatrix
     */
    Matrix<T> submat_del(const uint & p, const uint & q) const;

    /**
     * @brief compute determinant
     * matrix must be square
     * @return determinant
     */
    T det() const;

    /**
     * @brief computes minor of the matrix wrt row p and column q
     * (minor: determinant of the submatrix obtained by deletin
     * the p row and q column). The matrix must be square.
     * @param p row index
     * @param q column index
     * @return minor
    */
    T minor(const uint & p, const uint & q) const;

    /**
     * @brief computes cofactor of the matrix wrt row p and column q.
     * The matrix must be square.
     * @param p row index
     * @param q column index
     * @return cofactor
    */
    T cof(const uint & p, const uint & q) const;

    /**
     * @brief computes the cofactors matrix.
     * The input matrix must be square.
     * @return Matrix: cofactor matrix
    */
    Matrix<T> cof_mat() const;

    /**
     * @brief create adjoint of the matrix
     * 
     * @return Matrix: adjoint
     */
    Matrix<T> adj() const;

    /**
     * @brief compute inverse of the matrix using adjoint/det method
     * 
     * @return Matrix: inverse
     */
    Matrix<T> inv() const;

    /**
     * @brief compute left pseudo-inverse of the matrix
     * 
     * @return Matrix: left pseudo-inverse
     */
    Matrix<T> pinv_left() const;

    /**
     * @brief compute right pseudo-inverse of the matrix
     * 
     * @return Matrix: right pseudo-inverse
     */
    Matrix<T> pinv_right() const;
#pragma endregion matrix_operations

#pragma region decomposition_methods
    /**
     * @brief computes the reflector of the given vector
     * https://en.wikipedia.org/wiki/QR_decomposition
     * 
     * @return reflector vector
    */
    Matrix<T> reflector() const;

    /**
     * @brief Compute QR decomposition of the given matrix: A=Q*R 
     *        with Q orthogonal matrix and R upper triangular matrix
     *      http://matlab.izmiran.ru/help/techdoc/ref/mldivide.html
     *      https://rpubs.com/aaronsc32/qr-decomposition-householder
     * 
     * @param Q orthogonal matrix
     * @param R upper triangular matrix
     */
    void qr_dec(Matrix<T> & Q, Matrix<T> & R) const;

    /**
     * @brief Compute QR decomposition of the given matrix with partial permutation
     *      A*P=Q*R with: 
     *      Q orthogonal matrix 
     *      R upper triangular matrix
     *      P permutation matrix
     * 
     * @param Q orthogonal matrix
     * @param R upper triangular matrix
     * @param P permutation matrix
     */
    void qrp_dec(Matrix<T> & Q, Matrix<T> & R, Matrix<double> & P) const;


    /**
     * @brief Compute LUP decomposition of the given matrix: PA = LU with
     *        L lower triangular matrix and U upper triangular matrix and
     *        P permutation matrix
     *      https://www.tutorialspoint.com/cplusplus-program-to-perform-lu-decomposition-of-any-matrix
     * 
     * @param L lower triangular matrix
     * @param U upper triangular matrix
     * @param P permutation matrix
     * @return number of swap performed during permutation
     */
    uint lup_dec(Matrix<T> & L, Matrix<T> & U, Matrix<double> & P) const;

    /**
     * @brief preform hessenberg decomposition of the given matrix
     *        QHQ* = A
     *      https://en.wikipedia.org/wiki/Hessenberg_matrix
     * 
     * @param Q unitary matrix
     * @param H hessenberg matrix
    */
    void hessenberg_dec(Matrix<T> & Q, Matrix<T> & H) const;

#pragma endregion decomposition_methods

#pragma region eigen
    
    Matrix<T> implicit_double_QR_step() const;

    Matrix<c_double> get_eigenvalues(uint max_iterations = 1000, double tolerance = 1e-16) const;

    /**
     * @brief solution of the eigendecomposition porblem based on the double shift
     *      implicit QR algorithm with deflation.
     *      To assert the convergence, the algorithm check if during iteration
     *      the matrix is approaching an upper triangular matrix.
     * @param D matrix where the eigenvalues are saved as a vector
     * @param V matrix of the eigenvectors, saved as column vectors. The eigenvector on
     *      the i-th column is associated to the i-th eigenvalue of D
     * @param max_iterations maximum number of iterations of the algorithm
     * @param tolerance parameter to define the convergence criteria (elements of the matrix
     *      lower than tolerance are considered to be zero)
    */
    void eigen_dec(Matrix<c_double> & D, Matrix<c_double> & V, uint max_iterations = 1000, double tolerance = 1e-16) const;

#pragma endregion eigen

};

#pragma region ls_solution
    /**
     * @brief Solve the system U*x=B using a backward substitution algorithm.
     *        U must be an upper triangular square matrix (n*n) and B is the 
     *        known terms matrix with number of rows equalt to U -> B is (n*c_b)
     * 
     * @param U     upper triangular matrix (n*n)
     * @param B     known terms matrix (n*c_b)
     */
    template<typename T, typename V, typename W = RetType_t<T,V>>
    Matrix<W> backward_sub(Matrix<T> const & U, Matrix<V> const & B);

    /**
     * @brief Solve the system L*x=B using a forward substitution algorithm.
     *        L must be a lower triangular square matrix (n*n) and B is the 
     *        known terms matrix with number of rows equalt to L -> B is (n*c_b)
     * 
     * @param L     lower triangular matrix (n*n)
     * @param B     known terms matrix (n*c_b)
     */
    template<typename T, typename U, typename V = RetType_t<T,U>>
    Matrix<V> forward_sub(Matrix<T> const & L, Matrix<U> const & B);

    /**
     * @brief computes the left division A\B, which corresponds to solve the 
     *        linear equation system Ax=B or performing the operation (A^-1)*b.
     *        The rows of A must be equal the rows of B. The result has dimensions (r_a*c_b)
     *        The functions solves the problem depending on the dimensions of the 
     *        given matrices: 
     *        
     *        - A is a square matrix (r_a*c_r_common): A is decomposed using
     *          LUP decomposition: Ax=B -> PA = LU -> LUx = PB. Then, the problem is solved by
     *          subsequentially solving the two systems:
     *              - L*(U*x) = PB, thus solving it using forward substitution the system 
     *                L*y=PB equal to y=L\PB
     *              - U*x=y, thus solving it using backward substitution the system 
     *                U*x=y equal to x=U\y
     *        - A is rectangular and represent an overconstrained problem (rows > cols): A
     *          is decomposed using QRP decomposition (AP = QR) and performing the operation
     *          X = P*(R\(Q'*b)). Notably, using the householder projection for the QRP
     *          decomposition, matrix R is rectangular upper triangular, thus all rows below 
     *          the diagonal are zero. Such rows are discared along with the corresponding 
     *          rows of b. R\(Q'*b) is solved using backward substitution.
     * 
     * @param A             left hand division term
     * @param B             right hand division term
     * @return Matrix: result of the division
     * @throw invalid_argument if rows don't match
     */
    template<typename T, typename U, typename V = RetType_t<T,U>>
    Matrix<V> matrix_l_divide(Matrix<T> const & A, Matrix<U> const & B);

    /**
     * @brief solves the linear system A*x=B. 
     *        Calls matrix_l_divide.
     * @param A coefficient matrix (must be square nxn)
     * @param B known terms matrix (B.rows == A.cols)
     * @return Matrix: x (dimensions A.rows x B.cols)
     * @throw invalid_argument if rows don't match
     */
    template<typename T, typename U, typename V = RetType_t<T,U>>
    Matrix<V> solve_ls(Matrix<T> const & A, Matrix<U> const & B);

    /**
     * @brief computes the right division B/A translating it in a left 
     *        division problem following the equality B/A = (A.t\B.t).t
     *        where .t stands for transpose. 
     *        The columns of A must be equal the columns of B.
     *        The functions solves the problem depending on the dimensions of the 
     *        given matrices, as described in  matrix_l_divide.
     * 
     * @param B             left hand division term 
     * @param A             right hand division term
     * @return Matrix: result of the division
     * @throw invalid_argument if columns don't match
     */
    template<typename T, typename U, typename V = RetType_t<T,U>>
    Matrix<V> matrix_r_divide(Matrix<T> const & B, Matrix<U> const & A);

#pragma endregion ls_solution

#pragma region special_constructors
/**
 * @brief create identity matrix of shape dim*dim 
 * 
 * @param dim 
 * @return Matrix 
 */
Matrix<double> IdMat(const uint & dim);

/**
 * @brief create matrix of only ones of shape r*c
 * 
 * @param r rows
 * @param c columns
 * @return Matrix 
 */
Matrix<double> Ones(const uint & r, const uint & c);

/**
 * @brief create matrix with random values of shapre r*c
 * 
 * @param r rows
 * @param c columns
 * @return Matrix 
*/
template<typename T = double>
Matrix<T> RandMat(const uint & r, const uint & c);

/**
 * @brief given a vector of N elements, creates a NxN matrix
 *      with diagonal elements equal the vector elements.
 * @param v input vector 
*/
template<typename T>
Matrix<T> diag(const std::vector<T> & v);

/**
 * @brief given a row or column vector of N elements, creates a NxN matrix
 *      with diagonal elements equal the vector elements.
 * @param v input vector 
*/
template<typename T>
Matrix<T> diag(const Matrix<T>& v);
#pragma endregion special_constructors

std::ostream& operator<<(std::ostream& os, const uu_pair & p);

} // namespace MA

#include "linear_algebra_ma/matrices/matrices.tpp"
#include "linear_algebra_ma/matrices/operators.tpp"

#endif // MA_MATRICES_HPP
