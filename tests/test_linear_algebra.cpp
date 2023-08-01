#include <gtest/gtest.h>
#include "linear_algebra_ma/linearAlgebra.hpp"

using namespace MA;

using std::cout, std::endl, std::string;

TEST(Matrix, constructor_get) {
    EXPECT_NO_THROW(Matrix());
    Matrix m1, m2;
    EXPECT_EQ(m1.getR(), (uint)0);
    EXPECT_EQ(m1.getC(), (uint)0);
    EXPECT_EQ(m1.size(), (uint)0);
    EXPECT_EQ(m1.getV(), nullptr);

    EXPECT_NO_THROW(m1 = Matrix(2,2));
    EXPECT_EQ(m1.getR(), (uint)2);
    EXPECT_EQ(m1.getC(), (uint)2);
    EXPECT_EQ(m1.size(), (uint)4);
    EXPECT_NO_THROW(m1.getV()[3]);

    std::vector<double> v = {1,2,3,4,5,6};
    EXPECT_THROW(m1 = Matrix(4, 2, v), out_of_range);
    EXPECT_NO_THROW(m1 = Matrix(2, 3, v));
    EXPECT_EQ(m1.getR(), (uint)2);
    EXPECT_EQ(m1.getC(), (uint)3);
    EXPECT_EQ(m1.size(), (uint)6);
    EXPECT_EQ(m1(0,0), 1);
    EXPECT_EQ(m1(0,1), 2);
    EXPECT_EQ(m1(0,2), 3);
    EXPECT_EQ(m1(1,0), 4);
    EXPECT_EQ(m1(1,1), 5);
    EXPECT_EQ(m1(1,2), 6);

    EXPECT_NO_THROW(m2 = Matrix(m1));
    EXPECT_EQ(m1.getR(), (uint)2);
    EXPECT_EQ(m1.getC(), (uint)3);
    EXPECT_EQ(m1.size(), (uint)6);
    EXPECT_EQ(m2(0,0), 1);
    EXPECT_EQ(m2(0,1), 2);
    EXPECT_EQ(m2(0,2), 3);
    EXPECT_EQ(m2(1,0), 4);
    EXPECT_EQ(m2(1,1), 5);
    EXPECT_EQ(m2(1,2), 6);

    const double * tmp = m2.getV();
    EXPECT_EQ(tmp[0], 1);
    EXPECT_EQ(tmp[1], 2);
    EXPECT_EQ(tmp[2], 3);
    EXPECT_EQ(tmp[3], 4);
    EXPECT_EQ(tmp[4], 5);
    EXPECT_EQ(tmp[5], 6);
}

TEST(Matrix, access_operator){
    Matrix m2;
    Matrix m1(4,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25}
    );
    // single element
    EXPECT_THROW(m1(4,2), out_of_range);
    EXPECT_THROW(m1(1,5), out_of_range);
    EXPECT_EQ(m1(0,0), 10);
    EXPECT_EQ(m1(2,3), 21);

    // column vector
    EXPECT_THROW(m1({1,4},2), out_of_range);
    EXPECT_THROW(m1({1,3},5), out_of_range);
    EXPECT_THROW(m1({3,2},1), invalid_argument);
    EXPECT_NO_THROW(m2 = m1({1,3},1));
    EXPECT_EQ(m2.getR(), (uint)3);
    EXPECT_EQ(m2.getC(), (uint)1);
    EXPECT_EQ(m2(0,0), 15);
    EXPECT_EQ(m2(1,0), 19);
    EXPECT_EQ(m2(2,0), 23);
    EXPECT_NO_THROW(m2 = m1({2,2},1));
    EXPECT_EQ(m2.getR(), (uint)1);
    EXPECT_EQ(m2.getC(), (uint)1);
    EXPECT_EQ(m2(0,0), 19);

    // row vector
    EXPECT_THROW(m1(2,{1,4}), out_of_range);
    EXPECT_THROW(m1(5,{1,3}), out_of_range);
    EXPECT_THROW(m1(1,{3,2}), invalid_argument);
    EXPECT_NO_THROW(m2 = m1(1,{1,3}));
    EXPECT_EQ(m2.getR(), (uint)1);
    EXPECT_EQ(m2.getC(), (uint)3);
    EXPECT_EQ(m2(0,0), 15);
    EXPECT_EQ(m2(0,1), 16);
    EXPECT_EQ(m2(0,2), 17);
    EXPECT_NO_THROW(m2 = m1(1,{3,3}));
    EXPECT_EQ(m2.getR(), (uint)1);
    EXPECT_EQ(m2.getC(), (uint)1);
    EXPECT_EQ(m2(0,0), 17);

    // single index
    EXPECT_THROW(m1(4), out_of_range);
    EXPECT_EQ(m1(2), 20);
    m2 = m1({1,3},1);
    EXPECT_THROW(m2(3), out_of_range);
    EXPECT_EQ(m2(2), 23);
    m2 = m1(1,{1,3});
    EXPECT_THROW(m2(3), out_of_range);
    EXPECT_EQ(m2(2), 17);

    // double range access
    EXPECT_THROW(m1({1,4},{2,2}), out_of_range);
    EXPECT_THROW(m1({1,3},{1,5}), out_of_range);
    EXPECT_THROW(m1({4,1},{2,2}), invalid_argument);
    EXPECT_THROW(m1({1,3},{5,3}), invalid_argument);
    EXPECT_TRUE(m1({1,1},{3,3}) == Matrix(1,1,{17}));
    EXPECT_TRUE(m1({1,3},{3,3}) == Matrix(3,1,{17,21,25}));
    EXPECT_TRUE(m1({3,3},{1,2}) == Matrix(1,2,{23,24}));
    EXPECT_TRUE(m1({1,2},{1,2}) == Matrix(2,2,{15,16,19,20}));
    EXPECT_TRUE(m1({1,2},{0,3}) == Matrix(2,4,{14,15,16,17,18,19,20,21}));
}

TEST(Matrix, set) {
    Matrix m1(4,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25}
    );

    Matrix m2 = m1;
    std::vector<double> v = {0,1,2,3,
                             4,5,6,7,
                             8,9,10,11,
                             12,13,14,15};
    m2.setV(v);
    for(int i=0; i<4;i++) for(int j=0;j<4;j++) EXPECT_EQ(m2(i,j), v[j+i*4]);
    // cout << m2 << endl;

    m2 = m1;
    v = {0,1,2,
         3,4,5,
         6,7,8};
    EXPECT_THROW(m2.setV({1,4},{1,3}, v), out_of_range);
    EXPECT_THROW(m2.setV({1,2},{1,4}, v), out_of_range);
    EXPECT_THROW(m2.setV({2,1},{1,3}, v), invalid_argument);
    EXPECT_THROW(m2.setV({1,1},{3,1}, v), invalid_argument);
    EXPECT_THROW(m2.setV({0,3},{1,3}, v), out_of_range);
    EXPECT_NO_THROW(m2.setV({1,3},{1,3}, v));
    for(int i=0; i<4; i++) for(int j=0;j<4;j++){
        if(i==0 || j==0) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), v[j-1 + (i-1)*3]);
    }
    // cout << m2 << endl;

    m2 = m1;
    EXPECT_NO_THROW(m2.setV({1,2},{1,2}, v));
    for(int i=0; i<4; i++) for(int j=0;j<4;j++){
        if(i==0 || j==0 || i==3 || j==3) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), v[j-1 + (i-1)*2]);
    }
    // cout << m2 << endl;

    m2 = m1;
    EXPECT_NO_THROW(m2.setV({0,2},{2,3}, v));
    for(int i=0; i<4; i++) for(int j=0;j<4;j++){
        if(i>2 || j<2) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), v[j-2 + (i)*2]);
    }
    // cout << m2 << endl;

    Matrix m3(3,3,v);
    m2 = m1;
    EXPECT_THROW(m2.setV({1,4},{1,3}, m3), out_of_range);
    EXPECT_THROW(m2.setV({1,2},{1,4}, m3), out_of_range);
    EXPECT_THROW(m2.setV({2,1},{1,3}, m3), invalid_argument);
    EXPECT_THROW(m2.setV({1,1},{3,1}, m3), invalid_argument);
    EXPECT_THROW(m2.setV({0,3},{1,1}, m3), out_of_range);
    EXPECT_THROW(m2.setV({0,0},{0,4}, m3), out_of_range);
    EXPECT_NO_THROW(m2.setV({0,2},{2,3}, m3));
    for(int i=0; i<4; i++) for(int j=0;j<4;j++){
        if(i>2 || j<2) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), m3(i, j-2));
    }
    // cout << m2 << endl;
}

TEST(Matrix, transpose){
    Matrix m2;

    Matrix m1(3,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21}
    );
    EXPECT_NO_THROW(m2 = m1.t());
    EXPECT_EQ(m2.getR(), (uint)4);
    EXPECT_EQ(m2.getC(), (uint)3);
    for(int i=0; i<3; i++) for(int j=0;j<4;j++){
        EXPECT_EQ(m1(i,j), m2(j,i));
    }

    m1 = Matrix(1,5,
        {10,11,12,13,14}
    );
    EXPECT_NO_THROW(m2 = m1.t());
    EXPECT_EQ(m2.getR(), (uint)5);
    EXPECT_EQ(m2.getC(), (uint)1);
    for(int i=0; i<1; i++) for(int j=0;j<5;j++){
        EXPECT_EQ(m1(i,j), m2(j,i));
    }
}

TEST(Matrix, equal_operator){
    std::vector<double> v1 = 
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25};
    std::vector<double> v2 = 
        {10,11,12,13,
         4,5,6,7,
         8,9,10,11,
         12,13,14,15};

    EXPECT_FALSE(Matrix(3, 4, v1) == Matrix(4, 3, v1));
    EXPECT_FALSE(Matrix(3, 4, v1) == Matrix(3, 3, v1));
    EXPECT_TRUE(Matrix(2, 4, v1) == Matrix(2, 4, v1));
    EXPECT_TRUE(Matrix(2, 4, v1) == new Matrix(2, 4, v1));
    EXPECT_FALSE(Matrix(3, 4, v1) == Matrix(3, 4, v2));
    EXPECT_FALSE(Matrix(3, 4, v1) == new Matrix(3, 4, v2));
}

TEST(Matrix, unequal_operator){
    std::vector<double> v1 = 
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25};
    std::vector<double> v2 = 
        {10,11,12,13,
         4,5,6,7,
         8,9,10,11,
         12,13,14,15};

    EXPECT_TRUE(Matrix(3, 4, v1) != Matrix(4, 3, v1));
    EXPECT_TRUE(Matrix(3, 4, v1) != Matrix(3, 3, v1));
    EXPECT_FALSE(Matrix(2, 4, v1) != Matrix(2, 4, v1));
    EXPECT_FALSE(Matrix(2, 4, v1) != new Matrix(2, 4, v1));
    EXPECT_TRUE(Matrix(3, 4, v1) != Matrix(3, 4, v2));
    EXPECT_TRUE(Matrix(3, 4, v1) != new Matrix(3, 4, v2));
}

TEST(Matrix, sum_operator){
    std::vector<double> v1 =  {10,11,12,13,14,15,16,17,18,19,20,21};
    std::vector<double> v2 = v1;
    Matrix m1(3,4,v1);

    // int
    for(int i=0; i<v2.size(); i++) v2[i]+=7;
    EXPECT_TRUE(m1+7 == Matrix(3,4, v2));
    m1+=7;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)+=3);
    EXPECT_NO_THROW(Matrix(3,0)+3);

    // double
    m1 = Matrix(3,4,v1);
    v2 = v1;
    for(int i=0; i<v2.size(); i++) v2[i]+=(-3.14);
    EXPECT_TRUE(m1+(-3.14) == Matrix(3,4, v2));
    m1+=(-3.14);
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)+=7.72);
    EXPECT_NO_THROW(Matrix(3,0)+8.13);

    // matrix
    m1 = Matrix(3,4,v1);
    v2 = v1;
    Matrix m2 = m1 + 2.28;
    for(int i=0; i<v2.size(); i++) v2[i] = v2[i] + v2[i] + 2.28;
    EXPECT_TRUE((m1+m2) == Matrix(3,4, v2));
    m1+=m2;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));

    Matrix m3(4,3, v2);
    EXPECT_THROW(m1+m3, invalid_argument);
    EXPECT_THROW(m1+=m3, invalid_argument);
}

TEST(Matrix, subtract_operator){
    std::vector<double> v1 =  {10,11,12,13,14,15,16,17,18,19,20,21};
    std::vector<double> v2 = v1;
    Matrix m1(3,4,v1);

    // int
    for(int i=0; i<v2.size(); i++) v2[i]-=7;
    EXPECT_TRUE(m1-7 == Matrix(3,4, v2));
    m1-=7;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)-=3);
    EXPECT_NO_THROW(Matrix(3,0)-3);

    // double
    m1 = Matrix(3,4,v1);
    v2 = v1;
    for(int i=0; i<v2.size(); i++) v2[i]-=(-3.14);
    EXPECT_TRUE(m1-(-3.14) == Matrix(3,4, v2));
    m1-=(-3.14);
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)-=7.72);
    EXPECT_NO_THROW(Matrix(3,0)-8.13);

    // matrix
    m1 = Matrix(3,4,v1);
    v2 = v1;
    Matrix m2 = m1 + 2;
    for(int i=0; i<v2.size(); i++) v2[i] = v2[i] - v2[i] - 2;
    EXPECT_TRUE((m1-m2) == Matrix(3,4, v2));
    // cout << m1-m2 << endl;
    // cout << Matrix(3,4, v2) << endl;
    m1-=m2;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));

    Matrix m3(4,3, v2);
    EXPECT_THROW(m1-m3, invalid_argument);
    EXPECT_THROW(m1-=m3, invalid_argument);

}

TEST(Matrix, multiply_operator){
    std::vector<double> v1 =  {1,3,5,9,1,3,1,7,4,3,9,7};
    std::vector<double> v2 = v1;
    Matrix m1(3,4,v1);

    // int
    for(int i=0; i<v2.size(); i++) v2[i]*=2;
    EXPECT_TRUE(m1*2 == Matrix(3,4, v2));
    m1*=2;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)*=2);
    EXPECT_NO_THROW(Matrix(3,0)*2);

    // double
    m1 = Matrix(3,4,v1);
    v2 = v1;
    for(int i=0; i<v2.size(); i++) v2[i]*=2.5;
    EXPECT_TRUE(m1*2.5 == Matrix(3,4, v2));
    m1*=2.5;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)*=2.5);
    EXPECT_NO_THROW(Matrix(3,0)*2.5);

    // matrix
    m1 = Matrix(3,4,v1);
    v2 = v1;
    Matrix m2 = m1({1,2}, {0,3}).t();
    EXPECT_TRUE((m1*m2) == Matrix(3,2, {78,121,60,71,71,155}));
    m1*=m2;
    EXPECT_TRUE(m1 == Matrix(3,2, {78,121,60,71,71,155}));
    m1 = Matrix(3,4,v1);
    EXPECT_TRUE(m1(1,{1,3})*m1({0,2},2) == Matrix(1,1, {79}));
    m1*=m1.t();
    EXPECT_TRUE(m1 == Matrix(3,3, {116,78,121,78,60,71,121,71,155}));
    m1 = Matrix(3,4,v1);
    m2 = m1({1,2}, {0,3});
    EXPECT_THROW(m1*m2, invalid_argument);
    EXPECT_THROW(m1*=m1, invalid_argument);
}

TEST(Matrix, divide_operator){
    std::vector<double> v1 =  {1,3,5,9,1,3,1,7,4,3,9,7};
    std::vector<double> v2 = v1;
    Matrix m1(3,4,v1);

    // int
    for(int i=0; i<v2.size(); i++) v2[i]/=2;
    EXPECT_TRUE(m1/2 == Matrix(3,4, v2));
    m1/=2;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)/=2);
    EXPECT_NO_THROW(Matrix(3,0)/2);

    // double
    m1 = Matrix(3,4,v1);
    v2 = v1;
    for(int i=0; i<v2.size(); i++) v2[i]/=2.5;
    EXPECT_TRUE(m1/2.5 == Matrix(3,4, v2));
    m1/=2.5;
    EXPECT_TRUE(m1 == Matrix(3,4, v2));
    EXPECT_NO_THROW(Matrix(0,3)/=2.5);
    EXPECT_NO_THROW(Matrix(3,0)/2.5);

    // matrix
    Matrix A(4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,2,1,-1}
    );
    Matrix b(1,4, {2,4,1,3});
    // A not square
    EXPECT_THROW( b / A({0,3},{1,3}), runtime_error);
    // shape doesn't match
    EXPECT_THROW((b(0,{1,3}) / A), invalid_argument);
    // underdetermined
    Matrix C(4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,0,0,0}
    );
    EXPECT_THROW(b / C, runtime_error);
    // solution
    EXPECT_TRUE(b / A.t() == Matrix(1,4, {1,1,0,-1}));

    // A not square
    EXPECT_THROW( b/=A({0,3},{1,3}), runtime_error);
    // shape doesn't match
    EXPECT_THROW((b(0,{1,3})/=A), invalid_argument);
    // underdetermined
    C = Matrix(4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,0,0,0}
    );
    EXPECT_THROW(b/=C, runtime_error);
    // solution
    EXPECT_NO_THROW(b /= A.t());
    EXPECT_TRUE(b == Matrix(1,4, {1,1,0,-1}));

}

/*
TEST(Matrix, concatenate_operators){
    // TODO
}
*/

TEST(Matrix, submat_del){
    Matrix m2;

    Matrix m1(3,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21}
    );
    EXPECT_THROW(m1.submat_del(3,2), invalid_argument);
    EXPECT_THROW(m1.submat_del(0,4), invalid_argument);
    EXPECT_NO_THROW(m2 = m1.submat_del(1,2));
    EXPECT_EQ(m2.getR(), (uint)2);
    EXPECT_EQ(m2.getC(), (uint)3);
    for(int i=0; i<2; i++) for(int j=0;j<3;j++){
        uint k = i > 0 ? i+1 : i;
        uint l = j > 1 ? j+1 : j;
        EXPECT_EQ(m2(i,j), m1(k,l));
    }

    m1 = Matrix(1,5,
        {10,11,12,13,14}
    );
    EXPECT_NO_THROW(m2 = m1.submat_del(0,3));
    EXPECT_EQ(m2.getR(), (uint)0);
    EXPECT_EQ(m2.getC(), (uint)4);
    EXPECT_THROW(m2(0,0), out_of_range);

}

TEST(Matrix, determinant){
    Matrix m1(3,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21}
    );
    EXPECT_THROW(m1.det(), invalid_argument);

    m1 = Matrix(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    
    // 1x1
    EXPECT_EQ(m1({0,0},{0,0}).det(), 1);
    EXPECT_EQ(m1({2,2},{1,1}).det(), 3);
    // 2x2
    EXPECT_EQ(m1({1,2},{2,3}).det(), -56);
    // 3x3
    EXPECT_EQ(m1({1,3},{0,2}).det(), 110);
    // 4x4
    EXPECT_EQ(m1.det(), -376);
}

/*
TEST(Matrix, minor){
    // TODO
}

TEST(Matrix, cof){
    // TODO
}

TEST(Matrix, cof_mat){
    // TODO
}

TEST(Matrix, adj){
    // TODO
}

TEST(Matrix, inv){
    // TODO
}

TEST(Matrix, pinv_left){
    // TODO
}
*/

TEST(Matrix, norm2){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );

    EXPECT_THROW(m1.norm2(), invalid_argument);
    EXPECT_EQ(m1(2,{0,3}).norm2(), sqrt(155));
    EXPECT_EQ(m1({0,2},1).norm2(), sqrt(27));
    EXPECT_EQ(m1({3,3},{1,1}).norm2(), 2);
}

TEST(Matrix, normalize){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_THROW(m1.normalize(), invalid_argument);

    Matrix m2(1,4, {4,4,4,4});
    EXPECT_TRUE(m2.normalize() == Matrix(1,4, {0.5,0.5,0.5,0.5}));
}

TEST(Matrix, normalize_self){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_THROW(m1.normalize_self(), invalid_argument);

    Matrix m2(1,4, {4,4,4,4});
    m2.normalize_self();
    EXPECT_TRUE(m2 == Matrix(1,4, {0.5,0.5,0.5,0.5}));
}


/*
TEST(Matrix, qr_dec){
    // TODO
}
*/

TEST(Matrix, lu_dec){
    // matrix not square
    Matrix m1,L,U;
    m1 = Matrix (3,4,
        {1,3,5,9,
         0,2,1,7,
         4,1,8,2,
         5,2,1,9}
    );
    EXPECT_THROW(m1.lu_dec(L,U), invalid_argument);

    // matrix decomposable
    m1 = Matrix (4,4,
        {1,3,5,2,
         0,2,1,6,
         4,1,3,2,
         5,2,1,4}
    );
    EXPECT_NO_THROW(m1.lu_dec(L,U));
    EXPECT_TRUE(L*U == m1);

    // matrix not decomposable
    m1 = Matrix(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_THROW(m1.lu_dec(L,U), runtime_error);
}

TEST(Matrix, backward_sub){
    Matrix U,b;
    U = Matrix(4,4,
        {1,3,5,9,
         0,1,1,7,
         0,0,2,7,
         0,0,0,3}
    );
    b = Matrix(4,1, {2,6,1,3});

    // U not square
    EXPECT_THROW(Matrix::backward_sub(U({0,3},{1,3}), b), invalid_argument);
    // b.rows != U.cols
    EXPECT_THROW(Matrix::backward_sub(U, b({1,3},0)), invalid_argument);
    // underdetermined
    U(2,2) = 0;
    EXPECT_THROW(Matrix::backward_sub(U, b), runtime_error);
    U(2,2) = 2;
    
    // solve
    // cout << Matrix::backward_sub(U, b) << endl;
    EXPECT_TRUE(Matrix::backward_sub(U, b) == Matrix(4,1, {2,2,-3,1}));

}

TEST(Matrix, forward_sub){
    Matrix L,b;
    L = Matrix(4,4,
        {1,3,5,9,
         0,1,1,7,
         0,0,2,7,
         0,0,0,3}
    ).t();
    b = Matrix(4,1, {2,6,1,3});

    // L not square
    EXPECT_THROW(Matrix::forward_sub(L({0,3},{1,3}), b), invalid_argument);
    // b.rows != L.cols
    EXPECT_THROW(Matrix::forward_sub(L, b({1,3},0)), invalid_argument);
    // underdetermined
    L(2,2) = 0;
    EXPECT_THROW(Matrix::forward_sub(L, b), runtime_error);
    L(2,2) = 2;
    
    // solve
    // cout << Matrix::forward_sub(L, b) << endl;
    EXPECT_TRUE(Matrix::forward_sub(L, b) == Matrix(4,1, {2,0,-4.5,5.5}));

}

TEST(Matrix, matrix_l_divide){
    Matrix A,b,C;
    A = Matrix (4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,2,1,-1}
    );
    b = Matrix(4,1, {2,4,1,3});

    // A not square
    EXPECT_THROW(Matrix::matrix_l_divide(A({0,3},{1,3}), b), runtime_error);
    // shape doesn't match
    EXPECT_THROW(Matrix::matrix_l_divide(A, b({1,3},0)), invalid_argument);
    // underdetermined
    C = Matrix (4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,0,0,0}
    );
    EXPECT_THROW(Matrix::matrix_l_divide(C, b), runtime_error);

    // solution
    EXPECT_TRUE(Matrix::matrix_l_divide(A,b) == Matrix(4,1, {1,1,0,-1}));

    // Matrix L,U;
    // A.lu_dec(L,U);
    // cout << L << endl;
    // cout << U << endl;
    // cout << Matrix::matrix_l_divide(A,b) << endl;
}

TEST(Matrix, solve_ls){
    Matrix A,b,C;
    A = Matrix (4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,2,1,-1}
    );
    b = Matrix(4,1, {2,4,1,3});

    // A not square
    EXPECT_THROW(Matrix::solve_ls(A({0,3},{1,3}), b), invalid_argument);
    // shape doesn't match
    EXPECT_THROW(Matrix::solve_ls(A, b({1,3},0)), invalid_argument);
    // underdetermined
    C = Matrix (4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,0,0,0}
    );
    EXPECT_THROW(Matrix::solve_ls(C, b), runtime_error);

    // solution
    EXPECT_TRUE(Matrix::solve_ls(A,b) == Matrix(4,1, {1,1,0,-1}));
}

TEST(Matrix, matrix_r_divide){
    Matrix A,b,C;
    A = Matrix (4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,2,1,-1}
    );
    b = Matrix(1,4, {2,4,1,3});

    // A not square
    EXPECT_THROW(Matrix::matrix_r_divide(b, A({0,3},{1,3})), runtime_error);
    // shape doesn't match
    EXPECT_THROW(Matrix::matrix_r_divide(b(0,{1,3}), A), invalid_argument);
    // underdetermined
    C = Matrix (4,4,
        {1,3,4,2,
         0,2,1,-2,
         2,1,-3,2,
         0,0,0,0}
    );
    EXPECT_THROW(Matrix::matrix_r_divide(b, C), runtime_error);

    // solution
    EXPECT_TRUE(Matrix::matrix_r_divide(b,A.t()) == Matrix(1,4, {1,1,0,-1}));

    // Matrix L,U;
    // A.lu_dec(L,U);
    // cout << L << endl;
    // cout << U << endl;
    // cout << Matrix::matrix_r_divide(b,A.t()) << endl;
}