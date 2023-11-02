#ifndef MA_MATRICES
#define MA_MATRICES

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
#include <vector>
#include <array>

typedef unsigned int uint;

namespace MA
{

class Matrix {
protected:

    // number of decimals used in double comparison
    inline static uint double_precision = 10;
    // epsilon used in double comparison
    inline static double epsilon = std::pow(10,-10);

    uint _r;
    uint _c;
    double * _v;

public:

    /**
     * @brief set the double precision used in comparison
     * @param dp number of digits considered during comparison
     */
    inline static void set_double_precision(uint dp){ 
        if(dp > 15) std::cout << "WARNING: double precision very small, this may result in bad behavior" << std::endl;
        Matrix::double_precision = dp; 
        Matrix::epsilon = std::pow(10,-dp);
    }

    /**
     * @brief get the double precision used in comparison
     * @return number of digits considered during comparison
     */
    inline static uint get_double_precision(){ 
        return Matrix::double_precision; 
    }
    

#pragma region constructor_destructor
    Matrix();
    Matrix(const uint & r, const uint & c);
    Matrix(const uint & r, const uint & c, const std::vector<double> & v);
    Matrix(Matrix & m);
    Matrix(const Matrix & m);
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

    double& operator()(const uint & i);
    const double& operator()(const uint & i) const;
    Matrix operator()(std::pair<uint,uint> rs, const uint & c) const;
    Matrix operator()(const uint & r, std::pair<uint,uint> cs) const;
    Matrix operator()(std::pair<uint,uint> rs, std::pair<uint,uint> cs) const;
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix const * const m);

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
     * @brief uses the vecotr to fill in the elements of the matrix
     * @param v used vector
     * @throw out_of_range if v.size < this.size
     */
    void setV(const std::vector<double> & v);

    /**
     * @brief uses the vecotr to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns vector
     * @param v used vector
     * @throw out_of_range if v.size < submatrix size
     */
    void setV(std::pair<uint,uint> rs, std::pair<uint,uint> cs, const std::vector<double> & v);

    /**
     * @brief uses the given matrix to fill in the elements of submatrix
     * of this matrix defined by indeces rs and cs, extremes included
     * @param rs rows indeces
     * @param cs columns vector
     * @param m used matrix
     * @throw out_of_range if m has not enough rows or columns to fill in the 
     * submatrix
     */
    void setV(std::pair<uint,uint> rs, std::pair<uint,uint> cs, Matrix m);

    /**
     * @brief check if matrix is a vector
     * i.e. one of the dimension is 1
     * @return true if is at least one dim is == 1
     */
    bool is_vec() const;

    /**
     * @brief returns a matrix with the same underlaying
     * double array, but different shape. Size must be the same
     * @param r rows of new matrix
     * @param c cols of new matrix
     * @return Matrix: matrix with same data and specified shape
     */
    Matrix reshape(const uint & r, const uint & c) const;

    /**
     * @brief rehsape matrix into column vec
     * @return Matrix: column vec
     */
    Matrix to_c_vec() const;

    /**
     * @brief rehsape matrix into row vec
     * @return Matrix: row vec
     */
    Matrix to_r_vec() const;

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

#pragma region comparators
    bool operator==(const Matrix & m) const;
    bool operator==(const Matrix * const m) const;

    bool operator!=(const Matrix & m) const;
    bool operator!=(const Matrix * const m) const;
#pragma endregion comparators

#pragma region math_operators
    void operator+=(const double & k);
    void operator+=(const int & k);
    void operator+=(const Matrix & m);

    void operator-=(const double & k);
    void operator-=(const int & k);
    void operator-=(const Matrix & m);

    void operator*=(const double & k);
    void operator*=(const int & k);
    void operator*=(const Matrix & m);
    
    void operator/=(const double & k);
    void operator/=(const int & k);
    void operator/=(const Matrix & m);
#pragma endregion math_operators

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
    
    /**
     * @brief compute transpose
     * @return Matrix: transpose
     */
    Matrix t() const;

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
     * @brief test if the given matrix is singular
     * matrix must be square
     * @return bool: true if it's singular
     */
    bool is_sing() const;

    /**
     * @brief computes minor of the matrix wrt row p and column q
     * (minor: determinant of the submatrix obtained by deleting
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
     * using ((A'*A)^-1)*A'with A' transpose of A
     * 
     * @return Matrix: left pseudo-inverse
     */
    Matrix pinv_left() const;

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

#pragma region decomposition_methods
    /**
     * @brief Compute QR decomposition of the given matrix: A=Q*R 
     *        with Q orthogonal matrix and R upper triangular matrix
     *      //http://matlab.izmiran.ru/help/techdoc/ref/mldivide.html
     *      //https://rpubs.com/aaronsc32/qr-decomposition-householder
     * 
     * @param Q orthogonal matrix
     * @param R upper triangular matrix
     */
    void qr_dec(Matrix & Q, Matrix & R) const;


    //https://www.tutorialspoint.com/cplusplus-program-to-perform-lu-decomposition-of-any-matrix
    /**
     * @brief Compute LUP decomposition of the given matrix: PA = LU with
     *        L lower triangular matrix and U upper triangular matrix and
     *        P permutation matrix
     * 
     * @param L lower triangular matrix
     * @param U upper triangular matrix
     * @param P permutation matrix
     * @return number of swap performed during permutation
     */
    uint lup_dec(Matrix & L, Matrix & U, Matrix & P) const;
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
     *        linear equation system Ax=B. 
     *        The functions solves the problem depending on the dimensions of the 
     *        given matrices: 
     *        
     *        - A is a square matrix (r_a*c_r_common): The columns of A must be equal 
     *          the rows of B. The result has dimensions (r_a*c_b). A is decomposed using
     *          LUP decomposition: Ax=B -> PA = LU -> LUx = PB. Then, the problem is solved by
     *          subsequentially solving the two systems:
     *              - L*(U*x) = PB, thus solving it using forward substitution the system 
     *                L*y=PB equal to y=L\PB
     *              - U*x=y, thus solving it using backward substitution the system 
     *                U*x=y equal to x=U\y
     * 
     * @param A             left hand division term
     * @param B             right hand division term
     * @return Matrix: result of the division
     */
    static Matrix matrix_l_divide(Matrix const & A, Matrix const & B);

    /**
     * @brief solves the linear system A*x=B. 
     * Calls matrix_l_divide.
     * @param A coefficient matrix (must be square nxn)
     * @param B known terms matrix (B.rows == A.cols)
     * @return Matrix: x (dimensions A.rows x B.cols)
     */
    static Matrix solve_ls(Matrix const & A, Matrix const & B);

    /**
     * @brief computes the right division B/A translating it in a left 
     *        division problem following the equality B/A = (A.t\B.t).t
     *        where .t stands for transpose. 
     *        The functions solves the problem depending on the dimensions of the 
     *        given matrices, as described in  matrix_l_divide.
     * 
     * @param B             left hand division term 
     * @param A             right hand division term
     * @return Matrix: result of the division
     */
    static Matrix matrix_r_divide(Matrix const & B, Matrix const & A);

#pragma endregion ls_solution
};

// cout operators
std::ostream& operator<<(std::ostream& os, const Matrix& m);
std::ostream& operator<<(std::ostream& os, const Matrix * m);

#pragma region math_operators
Matrix operator+(const Matrix& m, const double& k);
Matrix operator+(const Matrix& m, const int & k);
Matrix operator+(const Matrix& m1, const Matrix& m2);

Matrix operator-(const Matrix& m, const double& k);
Matrix operator-(const Matrix& m, const int & k);
Matrix operator-(const Matrix& m1, const Matrix& m2);

Matrix operator*(const Matrix& m, const double& k);
Matrix operator*(const Matrix& m, const int & k);
Matrix operator*(const Matrix& m1, const Matrix& m2);

Matrix operator/(const Matrix& m, const double& k);
Matrix operator/(const Matrix& m, const int & k);
Matrix operator/(const Matrix& m1, const Matrix& m2);
#pragma endregion math_operators

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


// creation of particular matrices

/**
 * @brief create identity matrix of shape dim*dim 
 * 
 * @param dim 
 * @return Matrix 
 */
Matrix IdMat(const uint & dim);

/**
 * @brief create identity matrix of shape r*c
 * 
 * @param r rows 
 * @param c columns 
 * @return Matrix 
 */
Matrix IdMat(const uint & r, const uint & c);

/**
 * @brief create matrix of only ones of shape r*c
 * 
 * @param r 
 * @param c 
 * @return Matrix 
 */
Matrix Ones(const uint & r, const uint & c);

/**
 * @brief create matrix of only zeros of shape r*c
 * 
 * @param r 
 * @param c 
 * @return Matrix 
 */
Matrix Zeros(const uint & r, const uint & c);

/**
 * @brief create matrix of shape dim*dim with
 * diagonal elements equal to v
 * 
 * @param dim 
 * @return Matrix 
 */
Matrix diag(const uint & dim, const std::vector<double> & v);

/**
 * @brief create matrix of shape r*c with
 * diagonal elements equal to v
 * 
 * @param r 
 * @param c 
 * @return Matrix 
 */
Matrix diag(const uint & r, const uint & c, const std::vector<double> & v);

std::vector<double> diag(const Matrix & m);


} // namespace MA

#endif // MA_MATRICES