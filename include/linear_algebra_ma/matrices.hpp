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

using namespace std;
typedef unsigned int uint;

namespace MA
{

class Matrix {
protected: 
    uint _r;
    uint _c;
    double * _v;

public:

#pragma region constructor_destructor
    Matrix();
    Matrix(const uint & r, const uint & c);
    Matrix(const uint & r, const uint & c, std::vector<double> v);
    Matrix(Matrix & m);
    Matrix(const Matrix & m);
    ~Matrix();
#pragma endregion constructor_destructor

#pragma region get_set
    /**
     * @brief get operator using ()
     * 
     * @param r row index
     * @param c columns index
     * @return double& element in position (r,c)
     */
    double& operator()(const uint & r, const uint & c);

    /**
     * @brief const get operator using ()
     * 
     * @param r row index
     * @param c columns index
     * @return const double& element in position (r,c)
     */
    const double& operator()(const uint & r, const uint & c) const;

    double& operator()(const uint & i);
    const double& operator()(const uint & i) const;
    Matrix operator()(pair<uint,uint> rs, const uint & c) const;
    Matrix operator()(const uint & r, pair<uint,uint> cs) const;
    Matrix operator()(pair<uint,uint> rs, pair<uint,uint> cs) const;
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix const * const m);

    uint getR() const;
    uint getC() const;
    uint size() const;
    double const * getV() const;
    void setV(std::vector<double> v);
    void setV(pair<uint,uint> rs, pair<uint,uint> cs, std::vector<double> v);
    void setV(pair<uint,uint> rs, pair<uint,uint> cs, Matrix m);
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
     * 
     * @param m matrix to concatenate
     */
    void operator&=(const Matrix & m);

    /**
     * @brief concatenate matrices per rows
     * 
     * @param m matrix to concatenate
     */
    void operator|=(const Matrix & m);
    
    /**
     * @brief compute transpose
     * 
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
     * 
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
     * @brief Compute LU decomposition of the given matrix: A = L*U with
     *        L lower triangular matrix and U upper triangular matrix
     * 
     * @param L lower triangular matrix
     * @param U upper triangular matrix
     */
    void lu_dec(Matrix & L, Matrix & U) const;
#pragma endregion decomposition_methods

#pragma region ls_solution
    /**
     * @brief Solve the system U*x=B using a backward substitution algorithm.
     *        U must be an upper triangular square matrix (n*n) and B is the 
     *        known terms matrix with number of rows equalt to U -> B is (n*c_b)
     * 
     * @param U     upper triangular matrix (n*n)
     * @param B     known terms matrix (n*c_b)
     * @param res   matrix where to save the result (n*c_b)
     * @param n     dim of U
     * @param c_b   columns of B
     */
    static Matrix backward_sub(Matrix const & U, Matrix const & B);

    /**
     * @brief Solve the system L*x=B using a forward substitution algorithm.
     *        L must be a lower triangular square matrix (n*n) and B is the 
     *        known terms matrix with number of rows equalt to L -> B is (n*c_b)
     * 
     * @param L     lower triangular matrix (n*n)
     * @param B     known terms matrix (n*c_b)
     * @param res   matrix where to save the result (n*c_b)
     * @param n     dim of L
     * @param c_b   columns of B
     */
    static Matrix forward_sub(Matrix const & L, Matrix const & B);


    /**
     * @brief computes the left division A\B, which corresponds to solve the 
     *        linear equation system A*x=B. 
     *        The functions solves the problem depending on the dimensions of the 
     *        given matrices: 
     *        
     *        - A is a square matrix (r_a*c_r_common): The columns of A must be equal 
     *          the rows of B. The result has dimensions (r_a*c_b). A is decomposed using
     *          LU decomposition: A*x=B -> L*U*x = B. Then, the problem is solved by
     *          subsequentially solving the two systems:
     *              - L*(U*x) = B, thus solving it using forward substitution the system 
     *                L*par=B equal to par=L\B
     *              - U*x=par, thus solving it using backward substitution the system 
     *                U*x=par equal to x=U\par
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
ostream& operator<<(ostream& os, const Matrix& m);
ostream& operator<<(ostream& os, const Matrix * m);

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
 * @brief create matrix of only ones of shape r*c
 * 
 * @param r 
 * @param c 
 * @return Matrix 
 */
Matrix Ones(const uint & r, const uint & c);


Matrix diag(const uint & dim, double * v);
double * diag(const Matrix & m);


} // namespace MA

#endif // MA_MATRICES