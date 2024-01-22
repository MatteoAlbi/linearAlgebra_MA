#ifndef MA_MATRICES
#define MA_MATRICES

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cmath>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <random>

typedef unsigned int uint;
typedef std::pair<uint,uint> uu_pair;
#define ALL uu_pair{}

namespace MA
{

class Matrix {
protected:

    // number of decimals used in double comparison
    inline static uint double_precision = 10;
    // epsilon used in double comparison
    inline static double epsilon = 1e-10;
    // random number generator objects
    inline static std::uniform_real_distribution<double> unif{0.0, 1.0};
    inline static std::default_random_engine re;

    uint _r; // rows
    uint _c; // columns
    double * _v; // vector storing the elements per rows

public:

#pragma region static_methods
    /**
     * @brief set the double precision used in comparison
     * @param dp number of digits considered during comparison
     */
    inline static void set_double_precision(double dp = 10){ 
        if(dp > 15) std::cout << "WARNING: double precision very small, this may result in bad behavior" << std::endl;
        Matrix::double_precision = dp; 
        Matrix::epsilon = pow(10.0,-dp);
    }

    /**
     * @brief get the double precision used in comparison
     * @return margin of error considered during comparison
     */
    inline static double get_epsilon(){ 
        return Matrix::epsilon; 
    }

    /**
     * @brief generate a random number using the internally set seed
     * @return random number in range 0..1
     */
    inline static double rand(){
        return unif(re);
    }
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
    Matrix(const uint & r, const uint & c, std::vector<double> v);

    /**
     * @brief copy constructor
     * @param m matrix to copy
    */
    Matrix(Matrix & m);

    /**
     * @brief const copy constructor
     * @param m matrix to copy
    */
    Matrix(const Matrix & m);

    /**
     * @brief move constructor
    */
    Matrix(Matrix && m) noexcept;

    /**
     * @brief destructor, makes sure to free the memory
     *      allocated for v
    */
    ~Matrix();
#pragma endregion constructor_destructor

#pragma region get_set
    /**
     * @brief get operator using ()
     * @param r row index
     * @param c columns index
     * @return double& element in position (r,c)
     */
    double& operator()(const uint & r, const uint & c);

    /**
     * @brief const get operator using ()
     * @param r row index
     * @param c columns index
     * @return const double& element in position (r,c)
     */
    const double& operator()(const uint & r, const uint & c) const;

    /**
     * @brief if the object is a vector, access the i-th element,
     * if the object is a matrix, access the i-th element on the diagonal.
     * The returned value is a modifiable reference.
     * @param i index of access
     * @return element of the object
    */
    double& operator()(const uint & i);

    /**
     * @brief if the object is a vector, access the i-th element,
     * if the object is a matrix, access the i-th element on the diagonal.
     * The returned value is a constant refrence.
     * @param i index of access
     * @return element of the object
    */
    const double& operator()(const uint & i) const;

    /**
     * @brief extract the elements on rows rs and column c.
     * @param rs rows index pair (both extremes included)
     * @param c column index
     * @return submatrix of rows rs and column c
    */
    Matrix operator()(uu_pair rs, const uint & c) const;

    /**
     * @brief extract the elements on row r and columns cs.
     * @param r row index 
     * @param cs columns index pair (both extremes included)
     * @return submatrix of rows r and column cs
    */
    Matrix operator()(const uint & r, uu_pair cs) const;

    /**
     * @brief extract the elements on rows rs and columns cs.
     * @param rs rows index pair (both extremes included)
     * @param cs columns index pair (both extremes included)
     * @return submatrix of rows rs and columns cs
    */
    Matrix operator()(uu_pair rs, uu_pair cs) const;

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
    double const * v() const;

    /**
     * @brief uses the vector to fill in the elements of the matrix
     * @param v used vector
     * @throw out_of_range if v.size < this.size
     */
    void set(std::vector<double> v);

    /**
     * @brief uses the vector to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns vector
     * @param v used vector
     * @throw out_of_range if v.size < submatrix size
     */
    void set(uu_pair rs, uu_pair cs, std::vector<double> v);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs
     * @param r row index
     * @param c column index
     * @param m used matrix
     * @throw out_of_range if r or c are out of range
     */
    void set(uint r, uint c, double x);

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
    void set(uu_pair rs, uint c, Matrix m);

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
    void set(uint r, uu_pair cs, Matrix m);

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
    void set(uu_pair rs, uu_pair cs, Matrix m);

    /**
     * @brief returns a matrix with the same underlaying
     * double array, but different shape. Size must be the same
     * @param r rows of new matrix
     * @param c cols of new matrix
     * @return Matrix: matrix with same data and specified shape
     */
    Matrix reshape(const uint & r, const uint & c) const;

    /**
     * @brief extract diagonal of the matrix
     * @return vector (matrix with dim n,1) containing the diagonal elements
    */
    Matrix diag() const;

    /**
     * @brief swap the given rows of the matrix
     * @param r1 first row to swap
     * @param r2 second row to swap
     * @throw invalid argument if r1 or r2 exceeds the matric row indeces
     */
    void swap_rows(const uint & r1, const uint & r2);

    /**
     * @brief swap the given columns of the matrix
     * @param c1 first column to swap
     * @param c2 second column to swap
     * @throw invalid argument if c1 or c2 exceeds the matric column indeces
     */
    void swap_cols(const uint & c1, const uint & c2);
#pragma endregion get_set

#pragma region operators
    /**
     * @brief assignment operator from reference
     * @param m matrix to assign from
    */
    Matrix& operator=(const Matrix& m);

    /**
     * @brief assignment operator from pointer
     * @param m pointer to matrix to assign from
    */
    Matrix& operator=(Matrix const * const m);

    /**
     * @brief move assignment operator
     * @param m matrix to move
    */
    Matrix& operator=(Matrix && m) noexcept;

    bool operator==(const Matrix & m) const;
    bool operator==(const Matrix * const m) const;

    bool operator!=(const Matrix & m) const;
    bool operator!=(const Matrix * const m) const;

    void operator+=(const double & k);
    void operator+=(const Matrix & m);

    void operator-=(const double & k);
    void operator-=(const Matrix & m);

    void operator*=(const double & k);
    void operator*=(const Matrix & m);
    
    void operator/=(const double & k);
    void operator/=(const Matrix & m);

    /**
     * @brief concatenate matrices per columns
     * @param m matrix to concatenate
     */
    void operator&=(const Matrix & m);

    /**
     * @brief concatenate matrices per rows
     * @param m matrix to concatenate
     */
    void operator|=(const Matrix & m);
#pragma endregion operators

#pragma region vector
    /**
     * @brief reshape matrix into column vec
     * @return Matrix: column vec
     */
    Matrix to_c_vec() const;

    /**
     * @brief reshape matrix into row vec
     * @return Matrix: row vec
     */
    Matrix to_r_vec() const;

    /**
     * @brief computes dot product this.v
     * @param v second vector
     * @return double: dot product 
     */
    double dot(const Matrix & v) const;

    /**
     * @brief computes cross product this*v
     * implemented only for vectors of length 2 and 3
     * @param v second vector
     * @return Matrix: cross product (as column vec)
     */
    Matrix cross(const Matrix & v) const;

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
    Matrix normalize() const;

    /**
     * @brief normalized the vector, modifying current object
     * 
     */
    void normalize_self();
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
     * @brief compute transpose
     * @return Matrix: transpose
     */
    Matrix t() const;

    /**
     * @brief compute submatrix matrix obtained by deleting
     * the p-th row and q-t column
     * 
     * @param p row index
     * @param q column index
     * @return Matrix: submatrix
     */
    Matrix submat_del(const uint & p, const uint & q) const;

    /**
     * @brief compute determinant
     * matrix must be square
     * @return double: determinant
     */
    double det() const;

    /**
     * @brief computes minor of the matrix wrt row p and column q
     * (minor: determinant of the submatrix obtained by deletin
     * the p row and q column). The matrix must be square.
     * @param p row index
     * @param q column index
     * @return double: minor
    */
    double minor(const uint & p, const uint & q) const;

    /**
     * @brief computes cofactor of the matrix wrt row p and column q.
     * The matrix must be square.
     * @param p row index
     * @param q column index
     * @return double: cofactor
    */
    double cof(const uint & p, const uint & q) const;

    /**
     * @brief computes the cofactors matrix.
     * The input matrix must be square.
     * @return Matrix: cofactor matrix
    */
    Matrix cof_mat() const;

    /**
     * @brief create adjoint of the matrix
     * 
     * @return Matrix: adjoint
     */
    Matrix adj() const;

    /**
     * @brief compute inverse of the matrix using adjoint/det method
     * 
     * @return Matrix: inerse
     */
    Matrix inv() const;

    /**
     * @brief compute left pseudo-inverse of the matrix
     * 
     * @return Matrix: left pseudo-inverse
     */
    Matrix pinv_left() const;

    /**
     * @brief compute right pseudo-inverse of the matrix
     * 
     * @return Matrix: right pseudo-inverse
     */
    Matrix pinv_right() const;
#pragma endregion matrix_operations

#pragma region decomposition_methods
    /**
     * @brief Compute QR decomposition of the given matrix: A=Q*R 
     *        with Q orthogonal matrix and R upper triangular matrix
     *      http://matlab.izmiran.ru/help/techdoc/ref/mldivide.html
     *      https://rpubs.com/aaronsc32/qr-decomposition-householder
     * 
     * @param Q orthogonal matrix
     * @param R upper triangular matrix
     */
    void qr_dec(Matrix & Q, Matrix & R) const;

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
    void qrp_dec(Matrix & Q, Matrix & R, Matrix & P) const;


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
    uint lup_dec(Matrix & L, Matrix & U, Matrix & P) const;

    /**
     * @brief preform hessenberg decomposition of the given matrix
     *        QHQ* = A
     *      https://en.wikipedia.org/wiki/Hessenberg_matrix
     * 
     * @param Q unitary matrix
     * @param H hessenberg matrix
    */
    void hessenberg_dec(Matrix & Q, Matrix & H) const;

#pragma endregion decomposition_methods

#pragma region ls_solution
    /**
     * @brief Solve the system U*x=B using a backward substitution algorithm.
     *        U must be an upper triangular square matrix (n*n) and B is the 
     *        known terms matrix with number of rows equalt to U -> B is (n*c_b)
     * 
     * @param U     upper triangular matrix (n*n)
     * @param B     known terms matrix (n*c_b)
     */
    static Matrix backward_sub(Matrix const & U, Matrix const & B);

    /**
     * @brief Solve the system L*x=B using a forward substitution algorithm.
     *        L must be a lower triangular square matrix (n*n) and B is the 
     *        known terms matrix with number of rows equalt to L -> B is (n*c_b)
     * 
     * @param L     lower triangular matrix (n*n)
     * @param B     known terms matrix (n*c_b)
     */
    static Matrix forward_sub(Matrix const & L, Matrix const & B);


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
     *          X = P*(R\(Q'*b)). Notably, using the householder projectiob for the QRP
     *          decomposition, matrix R is rectangular upper triangular, thus all rows below 
     *          the diagonal are zero. Such rows are discared along with the corresponding 
     *          rows of b. R\(Q'*b) is solved using backward substitution.
     * 
     * @param A             left hand division term
     * @param B             right hand division term
     * @return Matrix: result of the division
     * @throw invalid_argument if rows don't match
     */
    static Matrix matrix_l_divide(Matrix const & A, Matrix const & B);

    /**
     * @brief solves the linear system A*x=B. 
     *        Calls matrix_l_divide.
     * @param A coefficient matrix (must be square nxn)
     * @param B known terms matrix (B.rows == A.cols)
     * @return Matrix: x (dimensions A.rows x B.cols)
     * @throw invalid_argument if rows don't match
     */
    static Matrix solve_ls(Matrix const & A, Matrix const & B);

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
    static Matrix matrix_r_divide(Matrix const & B, Matrix const & A);

#pragma endregion ls_solution

#pragma region eigen

    /**
     * @brief solution of the eigendecomposition porblem based on the QR algorithm.
     *      To assert the convergence, the algorithm check if during iteration
     *      the matrix is approaching an upper triangular matrix.
     * @param D matrix where the eigenvalues are saved as a vector
     * @param V matrix of the eigenvectors, saved as column vectors. The eigenvector on
     *      the i-th column is associated to the i-th eigenvalue of D
     * @param max_iterations maximum number of iterations of the algorithm
     * @param tolerance parameter to define the convergence criteria (elements of the matrix
     *      lower than tolerance are considered to be zero)
    */
    void eigen_QR(Matrix & D, Matrix & V, uint max_iterations = 1000, double tolerance = 1e-10) const;
    
    /**
     * @brief solution of the eigendecomposition porblem based on the QR algorithm with shift.
     *      To assert the convergence, the algorithm check if during iteration
     *      the matrix is approaching an upper triangular matrix.
     * @param D matrix where the eigenvalues are saved as a vector
     * @param V matrix of the eigenvectors, saved as column vectors. The eigenvector on
     *      the i-th column is associated to the i-th eigenvalue of D
     * @param max_iterations maximum number of iterations of the algorithm
     * @param tolerance parameter to define the convergence criteria (elements of the matrix
     *      lower than tolerance are considered to be zero)
    */
    void eigen_QR_shift(Matrix & D, Matrix & V, uint max_iterations = 1000, double tolerance = 1e-10) const;

#pragma endregion eigen

};

#pragma region operators
std::ostream& operator<<(std::ostream& os, const Matrix& m);
std::ostream& operator<<(std::ostream& os, const Matrix * m);
std::ostream& operator<<(std::ostream& os, const uu_pair & p);
std::ostream& operator<<(std::ostream& os, const uu_pair * p);

Matrix operator+(const Matrix& m, const double& k);
Matrix operator+(const double& k, const Matrix& m);
Matrix operator+(const Matrix& m1, const Matrix& m2);

Matrix operator-(const Matrix& m);
Matrix operator-(const Matrix& m, const double& k);
Matrix operator-(const double& k, const Matrix& m);
Matrix operator-(const Matrix& m1, const Matrix& m2);

Matrix operator*(const Matrix& m, const double& k);
Matrix operator*(const double& k, const Matrix& m);
Matrix operator*(const Matrix& m1, const Matrix& m2);

Matrix operator/(const Matrix& m, const double& k);
Matrix operator/(const double& k, const Matrix& m);
Matrix operator/(const Matrix& m1, const Matrix& m2);

/**
 * @brief concatenates matrices per columns
 * 
 * @param m1 first matrix
 * @param m2 second matrix
 * @return Matrix 
 */
Matrix operator&(const Matrix& m1, const Matrix& m2);

/**
 * @brief concatenates matrices per rows
 * 
 * @param m1 first matrix
 * @param m2 second matrix
 * @return Matrix
 */
Matrix operator|(const Matrix& m1, const Matrix& m2);
#pragma endregion operators

#pragma region special_constructors
/**
 * @brief create identity matrix of shape dim*dim 
 * 
 * @param dim 
 * @return Matrix 
 */
Matrix IdMat(const uint & dim);

/**
 * @brief create matrix of only ones of shape r*c
 * 
 * @param r rows
 * @param c columns
 * @return Matrix 
 */
Matrix Ones(const uint & r, const uint & c);

/**
 * @brief create matrix with random values of shapre r*c
 * 
 * @param r rows
 * @param c columns
 * @return Matrix 
*/
Matrix RandMat(const uint & r, const uint & c);

/**
 * @brief given a vector of N elements, creates a NxN matrix
 *      with diagonal elements equal the vector elements.
 * @param v input vector 
*/
Matrix diag(std::vector<double> v);

/**
 * @brief given a row or column vector of N elements, creates a NxN matrix
 *      with diagonal elements equal the vector elements.
 * @param v input vector 
*/
Matrix diag(const Matrix& v);
#pragma endregion special_constructors

} // namespace MA

#endif // MA_MATRICES