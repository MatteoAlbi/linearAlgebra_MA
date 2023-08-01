#include "linear_algebra_ma/linearAlgebra.hpp"

namespace MA
{

Matrix::Matrix(){
    this->_r = 0;
    this->_c = 0;
    this->_v = nullptr;
}

Matrix::Matrix(const uint & r, const uint & c) {
    this->_r = r;
    this->_c = c;
    this->_v = new double[this->size()]();
}

Matrix::Matrix(const uint & r, const uint & c, std::vector<double> v) {
    this->_r = r;
    this->_c = c;
    this->_v = new double[this->size()];
    this->setV(v);   
}

Matrix::Matrix(Matrix & m){
    this->_r = m._r;
    this->_c = m._c;
    this->_v = new double[m.size()];
    std::copy(m._v, m._v + m.size(), this->_v);
}

Matrix::Matrix(const Matrix & m){
    this->_r = m._r;
    this->_c = m._c;
    this->_v = new double[m.size()];
    std::copy(m._v, m._v + m.size(), this->_v);
}

Matrix::~Matrix() {
    if(this->_v != nullptr) delete[] this->_v;
}


uint Matrix::getR() const{return this->_r;}
uint Matrix::getC() const{return this->_c;}
uint Matrix::size() const{return this->_r * this->_c;}
double const * Matrix::getV() const{return this->_v;}


void Matrix::setV(std::vector<double> v){
    if(v.size() < this->size()) throw out_of_range("Not enough values in v to init the matrix");

    std::copy(v.begin(), v.begin() + this->size(), this->_v);
}

void Matrix::setV(pair<uint,uint> rs, pair<uint,uint> cs, std::vector<double> v){
    if(rs.second >= this->_r) throw out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw invalid_argument("Col first element must be <= of second");

    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(v.size() < nRows*nCols) throw out_of_range("Given vector doesn't have enough elements");

    for (uint i = 0; i < nRows; i++){
        for (uint j = 0; j < nCols; j++){
            this->operator()(i + rs.first, j + cs.first) = v[j + i * nCols];
        }
    }

}

void Matrix::setV(pair<uint,uint> rs, pair<uint,uint> cs, Matrix m){
    if(rs.second >= this->_r) throw out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw invalid_argument("Col first element must be <= of second");

    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(m.getR() < nRows) throw out_of_range("Given matrix doesn't have enough rows");
    if(m.getC() < nCols) throw out_of_range("Given matrix doesn't have enough cols");

    for (uint i = 0; i < nRows; i++){
        for (uint j = 0; j < nCols; j++){
            this->operator()(i + rs.first, j + cs.first) = m(i,j);
        }
    }

}


Matrix Matrix::t() const{
    Matrix ret = Matrix(this->_c, this->_r);
    for (uint i = 0; i<this->_r; i++) {
        for (uint j = 0; j<this->_c; j++) {
            ret._v[j * ret._c + i] = this->_v[i * this->_c + j];
        }
    }
    return ret;
}

Matrix Matrix::submat_del(const uint & p, const uint & q) const{
    if(this->_r <= p) throw invalid_argument("p greater than matrix rows");
    if(this->_c <= q) throw invalid_argument("q greater than matrix cols");

    uint i = 0, j = 0;

    Matrix ret = Matrix(this->_r-1,this->_c-1);

    // Looping for each element of the matrix
    for (uint row = 0; row < this->_r; row++){

        for (uint col = 0; col < this->_c; col++){

            // Copying into result matrix only those element
            // which are not in given row and column
            if (row != p && col != q){
                ret(i,j) = this->operator()(row,col);
                j++;
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == ret._c){
                    j = 0;
                    i++;
                }
            }
        }
    }

    return ret;
}

double Matrix::det() const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");

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

    else{
        double D = 0; // Initialize result

        Matrix sub;
        int sign = 1; // To store sign multiplier
    
        // Iterate for each element of first row
        for (uint j = 0; j < this->_c; j++){
            // submatrix of m->v[0][f]
            sub = this->submat_del(0, j);
            D += sign * this->_v[j] * sub.det();
    
            // terms are to be added with alternate sign
            sign = -sign;
        }

        return D;
    }
}

double Matrix::minor(const uint & p, const uint & q) const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");

    return this->submat_del(p,q).det();
}

double Matrix::cof(const uint & p, const uint & q) const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");

    int sign = (p+q % 2) == 0 ? 1 : -1;
    return this->submat_del(p,q).det() * sign;
}

Matrix Matrix::cof_mat() const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");

    Matrix ret(this->_r, this->_c);
    for(uint i=0; i<this->_r; i++) for(uint j=0; j<this->_c; j++){
        ret(i,j) = this->cof(i,j);
    }
}

Matrix Matrix::adj() const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");

    // same formula of cofactor matrix, but inverting indeces to get the transpose
    Matrix ret(this->_r, this->_c);
    for(uint i=0; i<this->_r; i++) for(uint j=0; j<this->_c; j++){
        ret(j,i) = this->cof(i,j);
    }
}

Matrix Matrix::inv() const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");

    // Find determinant of A
    double det = this->det(); 
    if (det == 0) throw runtime_error("Matrix not invertible: det = 0");
 
    // Find Inverse using formula "inverse(M) = adj(M)/det(M)"
    return this->adj() / det;
}

Matrix Matrix::pinv_left() const{
    try{
        return (this->t() * (*this)).inv() * this->t();
    }
    catch(const runtime_error& rte){
        std::cerr << rte.what() << '\n';
        throw runtime_error("Matrix (m.t * m) not invertible");
    }
    catch(const exception& e){
        std::cerr << e.what() << '\n';
        throw runtime_error("Unknown error occured");
    }
}

double Matrix::norm2() const{
    if(this->_r == 1 || this->_c == 1){
        double ret = 0;
        for(uint i=0; i< this->_r+this->_c-1; i++){
            ret += pow(this->_v[i],2);
        }
        return sqrt(ret);
    }
    else{
        // cout << "Norm only appliable to row/column vectors" << endl;
        // return NULL;
        throw invalid_argument("Norm only appliable to row/column vectors");
    }
}

Matrix Matrix::normalize() const{
    return Matrix(*this) / this->norm2();
}

void Matrix::normalize_self(){
    this->operator/=(this->norm2());
}



Matrix IdMat(const uint & dim){
    Matrix ret = Matrix(dim,dim);
    for(uint i=0; i<dim; i++){
        ret(i,i) = 1;
    }
    return ret;
}

Matrix Ones(const uint & r, const uint & c){
    return Matrix(r,c)+1;
}

Matrix diag(const uint & dim, double * v){
    Matrix ret = Matrix(dim,dim);
    
    for(uint i=0; i<dim; i++){
        ret(i,i) = v[i];
    }

    return ret;
}

double * diag(const Matrix & m){
    uint dim = min(m.getC(),m.getR());
    double * v = new double[dim]();

    for(uint i=0;i<dim; i++){
        v[i] = m(i,i);
    }

    return v;
}


#pragma region decomposition_methods

/*
void Matrix::qr_dec(Matrix & Q, Matrix & R) const{
    uint r = this->_r;
    uint c = this->_c;

    //set R to A
    R = this;
    Matrix Rtmp{r, c};
    
    uint n = std::min<double>(r, c);

    double v[r], tmp[r*r], hk[r*r];
    double H_list[n][r*r];
    double sign, vtv;
    uint v_dim;

    for(uint i=0; i<n; i++){
        //extract column of R
        Matrix v = R({i, r}, i);
        for(uint j=i; j<r; j++) v[j-i] = R[j*c + i];
        v_dim = r-i;

        //compute vk
        sign = v[0] < 0 ? -1 : 1;
        v[0] += sign * norm2(v, v_dim);

        normalize(v, v_dim);

        //compute H matrix
        matrix_mult_T(v, v, tmp, v_dim, v_dim, 1);
        for(uint j=0; j<v_dim*v_dim; j++) tmp[j] *= -2;
        for(uint j=0; j<v_dim; j++) tmp[j * v_dim + j] += 1;

        for(uint j=0; j<r; j++){
            for(uint k=0; k<r; k++){
                if(j >= i && k >= i) hk[j*r+k] = tmp[(j-i)*v_dim + (k-i)];
                else if(j == k) hk[j*r+k] = 1;
                else hk[j*r+k] = 0;
            }
        }

        //store H matrix to compute Q
        memcpy(H_list[i], hk, sizeof(double)*r*r);
        //update R
        matrix_mult(hk, R, Rtmp, r, c, r);
        memcpy(R, Rtmp, sizeof(double)*r*c);

    }

    //compute Q
    matrix_mult(H_list[0], H_list[1], Q, r, r, r);
    for(uint i=2; i<n; i++){
        matrix_mult(Q, H_list[i], tmp, r, r, r);
        memcpy(Q, tmp, sizeof(double)*r*r);
    }

}
*/

void Matrix::lu_dec(Matrix & L, Matrix & U) const{
    if(this->_r != this->_c) throw invalid_argument("The matrix must be square");
    uint n = this->_r;

    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++) {
            if (j < i) L[j*n + i] = 0;
            else {
                L[j*n + i] = A[j*n + i];
                for (uint k = 0; k < i; k++) L[j*n + i] -= L[j*n + k] * U[k*n + i];
            }
        }
        for (uint j = 0; j < n; j++) {
            if (j < i) U[i*n + j] = 0;
            else if (j == i) U[i*n + j] = 1;
            else {
                U[i*n + j] = A[i*n + j] / L[i*n + i];
                for (uint k = 0; k < i; k++)  U[i*n + j] -= L[i*n + k] * U[k*n + j] / L[i*n + i];
            }
        }
    }
}

#pragma endregion decomposition_methods


#pragma region ls_solution
// void backward_sub(const double *const U, const double *const B, 
//                  double *const res, const uint n, const uint c_b){
            
//     double tmp;
//     // printf("start i loop\n");
//     for(uint i=0; i<c_b; i++){ //col of res == col of B
//         // printf("start j loop with i:%d\n", i);
//         for(int j=n-1; j>=0; j--){ //row of res == row of U == row of B
//             tmp = 0;
//             // printf("start k loop with i:%d, j:%d\n", i, j);
//             for(int k=n-1; k>j; k--){ //col of U = row of res
//                 tmp += U[j*n + k] * res[k*c_b + i];
//                 // printf("i:%d j:%d k:%d\n",i,j,k);
//             }
            
//             res[j*c_b + i] = (B[j*c_b + i] - tmp) / U[j*n + j];
//             // printf("%f\t%f\n",tmp,res[j*c_b + i]);
//         }
//         // printf("\n");
//     }
    
// }


// void forward_sub(const double *const L, const double *const B, 
//                  double *const res, const uint n, const uint c_b){
//     double tmp;
    
//     for(uint j=0; j<n; j++){ //row of res == row of L == row of B
//         // printf("j: %d\n", j);
//         for(uint i=0; i<c_b; i++){ //col of res == col of B
//             // printf("i: %d\n", i);
//             tmp = 0;
//             for(uint k=0; k<j; k++){ //col of L = row of res
//                 // printf("k: %d\n", k);
//                 tmp += L[j*n + k] * res[k*c_b + i];
//             }
//             // printf("%f\n", tmp);
//             res[j*c_b + i] = (B[j*c_b + i] - tmp) / L[j*n + j];
//         }
//     }
    
// }


// void matrix_l_divide(const double *const A, const double *const B, 
//                  double *const res, const uint r_a, const uint c_b, const uint c_r_common){
    
//     if(r_a == c_r_common){ //square A
//         // printf("matrix A is square, using LU decomposition\n");
//         double L[r_a*r_a], U[r_a*r_a], res_tmp[r_a*c_b];

//         lu_dec(A, r_a, L, U);
//         forward_sub(L, B, res_tmp, r_a, c_b);
//         backward_sub(U, res_tmp, res, r_a, c_b);
//     }

// }


// void matrix_r_divide(const double *const B, const double *const A, 
//                  double *const res, const uint r_b, const uint c_a, const uint c_r_common){
    
//     double B_t[r_b*c_r_common], A_t[c_r_common*c_a], res_t[r_b*c_a];

//     transpose(B, B_t, r_b, c_r_common);
//     transpose(A, A_t, c_r_common, c_a);
    
//     matrix_l_divide(A_t, B_t, res_t, c_a, r_b, c_r_common);

//     transpose(res_t, res, c_a, r_b);
// }
#pragma endregion ls_solution

} // namespace MA

