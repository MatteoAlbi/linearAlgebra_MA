#include "linear_algebra_ma/matrices.hpp"

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


uint Matrix::r() const{return this->_r;}
uint Matrix::c() const{return this->_c;}
uint Matrix::size() const{return this->_r * this->_c;}
double const * Matrix::v() const{return this->_v;}


void Matrix::setV(std::vector<double> v){
    if(v.size() < this->size()) throw std::out_of_range("Not enough values in v to init the matrix");

    std::copy(v.begin(), v.begin() + this->size(), this->_v);
}

void Matrix::setV(std::pair<uint,uint> rs, std::pair<uint,uint> cs, std::vector<double> v){
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    if(rs.first == 0 && rs.second == 0) rs.second = this->_r-1;
    if(cs.first == 0 && cs.second == 0) cs.second = this->_c-1;

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(v.size() < nRows*nCols) throw std::out_of_range("Given vector doesn't have enough elements");

    for (uint i = 0; i < nRows; i++){
        for (uint j = 0; j < nCols; j++){
            this->operator()(i + rs.first, j + cs.first) = v[j + i * nCols];
        }
    }

}

void Matrix::setV(std::pair<uint,uint> rs, std::pair<uint,uint> cs, Matrix m){
    if(rs.second >= this->_r) throw std::out_of_range("Row index greater than this matrix rows");
    if(cs.second >= this->_c) throw std::out_of_range("Col index greater than this matrix cols");
    if(rs.first > rs.second) throw std::invalid_argument("Row first element must be <= of second"); 
    if(cs.first > cs.second) throw std::invalid_argument("Col first element must be <= of second");

    uint nRows = rs.second - rs.first + 1;
    uint nCols = cs.second - cs.first + 1;
    if(m.r() < nRows) throw std::out_of_range("Given matrix doesn't have enough rows");
    if(m.c() < nCols) throw std::out_of_range("Given matrix doesn't have enough cols");

    for (uint i = 0; i < nRows; i++){
        for (uint j = 0; j < nCols; j++){
            this->operator()(i + rs.first, j + cs.first) = m(i,j);
        }
    }

}


bool Matrix::is_vec() const{
    return this->_r == 1 || this->_c == 1;
}

Matrix Matrix::reshape(const uint & r, const uint & c) const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    if(r*c != this->size()) throw std::invalid_argument("New matrix size must match the current one");

    Matrix ret(r,c);
    std::copy(_v, _v + this->size(), ret._v);
    return ret;
}

Matrix Matrix::to_c_vec() const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    return this->reshape(this->size(), 1);
}

Matrix Matrix::to_r_vec() const{
    if(this->size() == 0) throw std::runtime_error("Matrix size is null");
    return this->reshape(1, this->size());
}

void Matrix::swap_rows(const uint & r1, const uint & r2){
    if(r1 >= _r || r2 >= _r) throw std::out_of_range("Given parameters exceed the matrix rows indeces");

    // extract r1
    Matrix tmp = this->operator()(r1, {0, _c-1});
    // substitute r2 into r1
    this->setV({r1,r1}, {0, _c-1}, this->operator()(r2, {0, _c-1}));
    // substitute r1 into r2
    this->setV({r2,r2}, {0, _c-1}, tmp);
}

void Matrix::swap_cols(const uint & c1, const uint & c2){
    if(c1 >= _c || c2 >= _c) throw std::out_of_range("Given parameters exceed the matrix columns indeces");

    // extract c1
    Matrix tmp = this->operator()({0, _r-1}, c1);
    // substitute c2 into c1
    this->setV({0, _r-1}, {c1,c1}, this->operator()({0, _r-1}, c2));
    // substitute c1 into c2
    this->setV({0, _r-1}, {c2,c2}, tmp);
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

double Matrix::dot(const Matrix & v) const{
    if(! this->is_vec()) throw std::invalid_argument("This obj is not a vector");
    if(! v.is_vec()) throw std::invalid_argument("Given obj is not a vector");
    if(this->size() != v.size()) throw std::invalid_argument("Vectors length don't match");

    double ret = 0;
    for(uint i=0; i<v.size(); i++) ret += this->operator()(i) * v(i);
    return ret;
}

Matrix Matrix::cross(const Matrix & v) const{
    if(! this->is_vec()) throw std::invalid_argument("This obj is not a vector");
    if(! v.is_vec()) throw std::invalid_argument("Given obj is not a vector");
    if(v.size() != 3 || this->size() != v.size()) throw std::invalid_argument("Cross product defined only for 3d vectors");

    Matrix ret(3,1);
    ret(0) =   ((*this)(1) * v(2)) - ((*this)(2) * v(1));
    ret(1) = -(((*this)(0) * v(2)) - ((*this)(2) * v(0)));
    ret(2) =   ((*this)(0) * v(1)) - ((*this)(1) * v(0));
    return ret;
}

Matrix Matrix::submat_del(const uint & p, const uint & q) const{
    if(this->_r <= p) throw std::invalid_argument("p greater than matrix rows");
    if(this->_c <= q) throw std::invalid_argument("q greater than matrix cols");

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
        Matrix L, U, P;
        uint n_swaps = this->lup_dec(L,U,P);
        double determinant = std::pow(-1.0, n_swaps);
        for(uint i=0; i<U.c(); i++) determinant *= U(i,i);
        return determinant;
    }
}

bool Matrix::is_sing() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");
    return this->det() == 0;
}

double Matrix::minor(const uint & p, const uint & q) const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    return this->submat_del(p,q).det();
}

double Matrix::cof(const uint & p, const uint & q) const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    int sign = (p+q % 2) == 0 ? 1 : -1;
    return this->submat_del(p,q).det() * sign;
}

Matrix Matrix::cof_mat() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    Matrix ret(this->_r, this->_c);
    for(uint i=0; i<this->_r; i++) for(uint j=0; j<this->_c; j++){
        ret(i,j) = this->cof(i,j);
    }

    return ret;
}

Matrix Matrix::adj() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    // same formula of cofactor matrix, but inverting indeces to get the transpose
    Matrix ret(this->_r, this->_c);
    for(uint i=0; i<this->_r; i++) for(uint j=0; j<this->_c; j++){
        ret(j,i) = this->cof(i,j);
    }

    return ret;
}

Matrix Matrix::inv() const{
    if(this->_r != this->_c) throw std::invalid_argument("The matrix must be square");

    return matrix_l_divide(*this, IdMat(_r));
}

Matrix Matrix::pinv_left() const{
    try{
        return (this->t() * (*this)).inv() * this->t();
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
        throw std::invalid_argument("Norm only appliable to row/column vectors");
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
    uint dim = std::min(m.c(),m.r());
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


uint Matrix::lup_dec(Matrix & L, Matrix & U, Matrix & P) const{
    if(_r != _c) throw std::invalid_argument("The matrix must be square");

    U = *this;
    L = IdMat(_r);
    P = IdMat(_r);
    uint ret = 0;

    for (uint i = 0; i < _c; i++) { // scroll columns
        // pivoting
        double u_max = 0;
        uint max_index = i;
        for (uint j = i; j < _r; j++){ // scroll rows
            // find max value
            if(u_max < std::fabs(U(j,i))){
                max_index = j;
                u_max = std::fabs(U(j,i));
            }
        }
        // if max value is zero we can skip this iteration
        if(u_max == 0) continue;

        // swap rows
        if(max_index != i){
            P.swap_rows(i, max_index);
            U.swap_rows(i, max_index);
            ret++;
        }

        // elimination
        for (uint j = i+1; j < _r; j++){ // scroll rows
            // compute l(j,i)
            U(j,i) = U(j,i) / U(i,i);
            // apply elimination over all elements of the row
            for (uint k = i+1; k < _c; k++){ // scroll columns
                U(j,k) -= U(j,i) * U(i,k);
            }
        }    
    }

    // compose matrices
    for(uint i=1; i<_r; i++) { //scroll rows
        for(uint j=0; j<i; j++){ //scroll columns
            L(i,j) = U(i,j);
            U(i,j) = 0;
        }
    }

    return ret;
}

#pragma endregion decomposition_methods


#pragma region ls_solution
Matrix Matrix::backward_sub(Matrix const & U, Matrix const & B){
    if(U.c() != U.r()) throw std::invalid_argument("Coefficient matrix U must be square");
    if(B.r() != U.c()) throw std::invalid_argument("Rows of B must be equal to the columns of U");
    // check all U diagonal elements are different from zero
    for(uint i=0; i<U.c(); i++) if(U(i,i) == 0){
        throw std::runtime_error("System of equation is underdetermined");
    }

    Matrix res(U.r(), B.c());

    double tmp;
    // printf("start i loop\n");
    for(uint i=0; i<B.c(); i++){ //col of res == col of B
        // printf("start j loop with i:%d\n", i);
        for(int j=U.r()-1; j>=0; j--){ //row of res == row of U == row of B
            tmp = 0;
            // printf("start k loop with i:%d, j:%d\n", i, j);
            for(int k=U.r()-1; k>j; k--){ //col of U = row of res
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


Matrix Matrix::forward_sub(Matrix const & L, Matrix const & B){
    if(L.c() != L.r()) throw std::invalid_argument("Coefficient matrix L must be square");
    if(B.r() != L.c()) throw std::invalid_argument("Rows of B must be equal to the columns of L");
    // check all L diagonal elements are different from zero
    for(uint i=0; i<L.c(); i++) if(L(i,i) == 0){
        throw std::runtime_error("System of equation is underdetermined");
    }

    Matrix res(L.r(), B.c());
    
    double tmp;
    for(uint j=0; j<L.r(); j++){ //row of res == row of L == row of B
        // printf("j: %d\n", j);
        for(uint i=0; i<B.c(); i++){ //col of res == col of B
            // printf("i: %d\n", i);
            tmp = 0;
            for(uint k=0; k<j; k++){ //col of L = row of res
                // printf("k: %d\n", k);
                tmp += L(j,k) * res(k,i);
            }
            // printf("%f\n", tmp);
            res(j,i) = (B(j,i) - tmp) / L(j,j);
        }
    }
    
    return res;
}


Matrix Matrix::matrix_l_divide(Matrix const & A, Matrix const & B){
    
    if(A.c() == A.r()){ //square A
        // check rows of B
        if(A.c() != B.r()){
            throw std::invalid_argument("Rows of B=(" + std::to_string(B.r()) + ") \
            must be equal the columns of A=" + std::to_string(A.c()));
        }
        // printf("matrix A is square, using LU decomposition\n");
        Matrix L, U, P, res_tmp;

        A.lup_dec(L, U, P);
        res_tmp = forward_sub(L, P*B);
        return backward_sub(U, res_tmp);
    }
    else{
        throw std::runtime_error("Methods for not square A are yet to be implemented");
    }

}

Matrix Matrix::solve_ls(Matrix const & A, Matrix const & B){
    if(A.c() != A.r()) throw std::invalid_argument("Coefficient matrix A must be square");
    return Matrix::matrix_l_divide(A,B);
}

Matrix Matrix::matrix_r_divide(Matrix const & B, Matrix const & A){
    return matrix_l_divide(A.t(), B.t()).t();
}
#pragma endregion ls_solution

} // namespace MA

