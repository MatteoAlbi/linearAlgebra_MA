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
    
}