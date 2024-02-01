#include <gtest/gtest.h>
#include "linear_algebra_ma/matrices.hpp"

using namespace MA;
using namespace std::complex_literals;

using std::cout, std::endl, std::string, std::complex;
using std::invalid_argument, std::runtime_error, std::out_of_range;


TEST(Matrix, constructor_getter) {
    std::vector<double> v{1,2,3,4,5,6};
    std::vector<complex<double>> v_comp{1.0+1i,2.0+1i,3.0+1i,4.0+1i,5.0+1i,6.0+1i};

    // empty constructor
    EXPECT_NO_THROW(Matrix());
    Matrix empty{};
    EXPECT_EQ(empty.r(), (uint)0);
    EXPECT_EQ(empty.c(), (uint)0);
    EXPECT_EQ(empty.v(), nullptr);

    // zeros initialization
    EXPECT_NO_THROW(Matrix zeros(2,2));
    Matrix zeros(2,2);
    EXPECT_EQ(zeros.r(), (uint)2);
    EXPECT_EQ(zeros.c(), (uint)2);
    EXPECT_NO_THROW((void) zeros.v()[3]);
    for(uint i=0; i<2; ++i) for(uint j=0; j<2; ++j) EXPECT_EQ(zeros(i,j), 0);

    // vector initialization
    EXPECT_THROW(Matrix from_vec(4, 2, v), out_of_range);
    EXPECT_NO_THROW(Matrix from_vec(2, 3, v));
    Matrix from_vec(2, 3, v);
    EXPECT_EQ(from_vec.r(), (uint)2);
    EXPECT_EQ(from_vec.c(), (uint)3);
    EXPECT_NO_THROW((void) from_vec.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(from_vec(i,j), k);

    // copy constructor
    EXPECT_NO_THROW(Matrix mat_copy(from_vec));
    Matrix mat_copy(from_vec);
    EXPECT_EQ(mat_copy.r(), (uint)2);
    EXPECT_EQ(mat_copy.c(), (uint)3);
    EXPECT_NO_THROW((void) mat_copy.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(mat_copy(i,j), k);

    // move constructor
    Matrix mat_moved(std::move(from_vec));
    EXPECT_EQ(from_vec.r(), (uint)0);
    EXPECT_EQ(from_vec.c(), (uint)0);
    EXPECT_EQ(from_vec.v(), nullptr);
    EXPECT_EQ(mat_moved.r(), (uint)2);
    EXPECT_EQ(mat_moved.c(), (uint)3);
    EXPECT_NO_THROW((void) mat_moved.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(mat_moved(i,j), k);


    // -- complex matrix -- //
    // empty constructor
    EXPECT_NO_THROW(Matrix<complex<double>>());
    Matrix<complex<double>> empty_comp{};
    EXPECT_EQ(empty_comp.r(), (uint)0);
    EXPECT_EQ(empty_comp.c(), (uint)0);
    EXPECT_EQ(empty_comp.v(), nullptr);

    // zeros initialization
    EXPECT_NO_THROW(Matrix<complex<double>> zeros_comp(2,2));
    Matrix<complex<double>> zeros_comp(2,2);
    EXPECT_EQ(zeros_comp.r(), (uint)2);
    EXPECT_EQ(zeros_comp.c(), (uint)2);
    EXPECT_NO_THROW((void) zeros_comp.v()[3]);
    for(uint i=0; i<2; ++i) for(uint j=0; j<2; ++j) EXPECT_EQ(zeros_comp(i,j), 0.0);

    // vector initialization
    EXPECT_THROW(Matrix<complex<double>> comp_from_vec(4, 2, v_comp), out_of_range);
    EXPECT_NO_THROW(Matrix<complex<double>> comp_from_vec(2, 3, v_comp));
    Matrix<complex<double>> comp_from_vec(2, 3, v_comp);
    EXPECT_EQ(comp_from_vec.r(), (uint)2);
    EXPECT_EQ(comp_from_vec.c(), (uint)3);
    EXPECT_NO_THROW((void) comp_from_vec.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(comp_from_vec(i,j), (double)k + 1i);

    // copy constructor
    EXPECT_NO_THROW(Matrix mat_comp_copy(comp_from_vec));
    Matrix mat_comp_copy(comp_from_vec);
    EXPECT_EQ(mat_comp_copy.r(), (uint)2);
    EXPECT_EQ(mat_comp_copy.c(), (uint)3);
    EXPECT_NO_THROW((void) mat_comp_copy.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(mat_comp_copy(i,j), (double)k + 1i);

    // move constructor
    Matrix mat_moved_comp(std::move(comp_from_vec));
    EXPECT_EQ(comp_from_vec.r(), (uint)0);
    EXPECT_EQ(comp_from_vec.c(), (uint)0);
    EXPECT_EQ(comp_from_vec.v(), nullptr);
    EXPECT_EQ(mat_moved_comp.r(), (uint)2);
    EXPECT_EQ(mat_moved_comp.c(), (uint)3);
    EXPECT_NO_THROW((void) mat_moved_comp.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(mat_moved_comp(i,j), (double)k + 1i);


    // -- from double to complex -- //
    // vector initialization
    EXPECT_THROW(Matrix<complex<double>> comp_from_vec_double(4, 2, v), out_of_range);
    EXPECT_NO_THROW(Matrix<complex<double>> comp_from_vec_double(2, 3, v));
    Matrix<complex<double>> comp_from_vec_double(2, 3, v);
    EXPECT_EQ(comp_from_vec_double.r(), (uint)2);
    EXPECT_EQ(comp_from_vec_double.c(), (uint)3);
    EXPECT_NO_THROW((void) comp_from_vec_double.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(comp_from_vec_double(i,j), (double)k);

    // copy constructor
    EXPECT_NO_THROW(Matrix<complex<double>> mat_comp_copy_double(mat_copy));
    Matrix<complex<double>> mat_comp_copy_double(mat_copy);
    EXPECT_EQ(mat_comp_copy_double.r(), (uint)2);
    EXPECT_EQ(mat_comp_copy_double.c(), (uint)3);
    EXPECT_NO_THROW((void) mat_comp_copy_double.v()[3]);
    for(uint i=0, k=1; i<2; ++i, ++k) for(uint j=0; j<2; ++j, ++k) EXPECT_EQ(mat_comp_copy_double(i,j), (double)k);


    // -- from complex to double -- //
    // EXPECT_ANY_THROW(Matrix<double> tmp(mat_comp_copy_double));
    // -> compilation error: cannot convert ‘const std::complex<double>’ to ‘double’ in assignment
}

TEST(Matrix, assignment_operator) {
    Matrix m1(3,2, {1,2,3,4,5,6});
    const Matrix m2(3,2, {7,8,9,10,11,12});
    Matrix m3, m4(4,4);

    // copy from ref
    EXPECT_NO_THROW(m3 = m1);
    EXPECT_EQ(m3, m1);
    // copy from const ref
    EXPECT_NO_THROW(m3 = m2);
    EXPECT_EQ(m3, m2);
    // move
    EXPECT_NO_THROW(m3 = std::move(m4));
    EXPECT_EQ(m3, Matrix(4,4));
    // matrix is empty cause is moved by the constructor
    // (swapped with default constructed matrix)
    EXPECT_EQ(m4, Matrix());

    
    // -- complex matrix -- //
    Matrix<complex<double>> m1c(3,2, {1.0+1i,2.0+1i,3.0+1i,4.0+1i,5.0+1i,6.0+1i});
    const Matrix<complex<double>> m2c(3,2, {7.0+1i,8.0+1i,9.0+1i,10.0+1i,11.0+1i,12.0+1i});
    Matrix<complex<double>> m3c, m4c(4,4);

    // copy from ref
    EXPECT_NO_THROW(m3c = m1c);
    EXPECT_EQ(m3c, m1c);
    // copy from const ref
    EXPECT_NO_THROW(m3c = m2c);
    EXPECT_EQ(m3c, m2c);
    // move
    EXPECT_NO_THROW(m3c = std::move(m4c));
    EXPECT_EQ(m3c, Matrix(4,4));
    // matrix is empty cause is moved by the constructor
    // (swapped with default constructed matrix)
    EXPECT_EQ(m4c, Matrix());


    // -- from double to complex -- //
    // copy from ref
    EXPECT_NO_THROW(m3c = m1);
    EXPECT_EQ(m3c, m1);
    // copy from const ref
    EXPECT_NO_THROW(m3c = m2);
    EXPECT_EQ(m3c, m2);
    // move
    EXPECT_NO_THROW(m3c = std::move(m3));
    // there is no move constructor for matrices with different typenames
    // -> m3 will not be modified
    EXPECT_EQ(m3c, m3);


    // -- from complex to double -- //
    // EXPECT_ANY_THROW(m1 = m1c);
    // -> compilation error: cannot convert ‘const std::complex<double>’ to ‘double’ in assignment
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
    EXPECT_EQ(m2.r(), (uint)3);
    EXPECT_EQ(m2.c(), (uint)1);
    EXPECT_EQ(m2(0,0), 15);
    EXPECT_EQ(m2(1,0), 19);
    EXPECT_EQ(m2(2,0), 23);
    EXPECT_NO_THROW(m2 = m1({2,2},1));
    EXPECT_EQ(m2.r(), (uint)1);
    EXPECT_EQ(m2.c(), (uint)1);
    EXPECT_EQ(m2(0,0), 19);

    // row vector
    EXPECT_THROW(m1(2,{1,4}), out_of_range);
    EXPECT_THROW(m1(5,{1,3}), out_of_range);
    EXPECT_THROW(m1(1,{3,2}), invalid_argument);
    EXPECT_NO_THROW(m2 = m1(1,{1,3}));
    EXPECT_EQ(m2.r(), (uint)1);
    EXPECT_EQ(m2.c(), (uint)3);
    EXPECT_EQ(m2(0,0), 15);
    EXPECT_EQ(m2(0,1), 16);
    EXPECT_EQ(m2(0,2), 17);
    EXPECT_NO_THROW(m2 = m1(1,{3,3}));
    EXPECT_EQ(m2.r(), (uint)1);
    EXPECT_EQ(m2.c(), (uint)1);
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
    EXPECT_EQ(m1({1,1},{3,3}), Matrix(1,1,{17}));
    EXPECT_EQ(m1({1,3},{3,3}), Matrix(3,1,{17,21,25}));
    EXPECT_EQ(m1({3,3},{1,2}), Matrix(1,2,{23,24}));
    EXPECT_EQ(m1({1,2},{1,2}), Matrix(2,2,{15,16,19,20}));
    EXPECT_EQ(m1({1,2},{0,3}), Matrix(2,4,{14,15,16,17,18,19,20,21}));


    // -- complex matrix -- //
    Matrix<complex<double>> m2c;
    Matrix<complex<double>> m1c(4,4,
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
         14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,
         22.0+1i,23.0+1i,24.0+1i,25.0+1i}
    );
    // single element
    EXPECT_THROW(m1c(4,2), out_of_range);
    EXPECT_THROW(m1c(1,5), out_of_range);
    EXPECT_EQ(m1c(0,0), 10.0+1i);
    EXPECT_EQ(m1c(2,3), 21.0+1i);

    // column vector
    EXPECT_THROW(m1c({1,4},2), out_of_range);
    EXPECT_THROW(m1c({1,3},5), out_of_range);
    EXPECT_THROW(m1c({3,2},1), invalid_argument);
    EXPECT_NO_THROW(m2c = m1c({1,3},1));
    EXPECT_EQ(m2c.r(), (uint)3);
    EXPECT_EQ(m2c.c(), (uint)1);
    EXPECT_EQ(m2c(0,0), 15.0+1i);
    EXPECT_EQ(m2c(1,0), 19.0+1i);
    EXPECT_EQ(m2c(2,0), 23.0+1i);
    EXPECT_NO_THROW(m2c = m1c({2,2},1));
    EXPECT_EQ(m2c.r(), (uint)1);
    EXPECT_EQ(m2c.c(), (uint)1);
    EXPECT_EQ(m2c(0,0), 19.0+1i);

    // row vector
    EXPECT_THROW(m1c(2,{1,4}), out_of_range);
    EXPECT_THROW(m1c(5,{1,3}), out_of_range);
    EXPECT_THROW(m1c(1,{3,2}), invalid_argument);
    EXPECT_NO_THROW(m2c = m1c(1,{1,3}));
    EXPECT_EQ(m2c.r(), (uint)1);
    EXPECT_EQ(m2c.c(), (uint)3);
    EXPECT_EQ(m2c(0,0), 15.0+1i);
    EXPECT_EQ(m2c(0,1), 16.0+1i);
    EXPECT_EQ(m2c(0,2), 17.0+1i);
    EXPECT_NO_THROW(m2c = m1c(1,{3,3}));
    EXPECT_EQ(m2c.r(), (uint)1);
    EXPECT_EQ(m2c.c(), (uint)1);
    EXPECT_EQ(m2c(0,0), 17.0+1i);

    // single index
    EXPECT_THROW(m1c(4), out_of_range);
    EXPECT_EQ(m1c(2), 20.0+1i);
    m2c = m1c({1,3},1);
    EXPECT_THROW(m2c(3), out_of_range);
    EXPECT_EQ(m2c(2), 23.0+1i);
    m2c = m1c(1,{1,3});
    EXPECT_THROW(m2c(3), out_of_range);
    EXPECT_EQ(m2c(2), 17.0+1i);

    // double range access
    EXPECT_THROW(m1c({1,4},{2,2}), out_of_range);
    EXPECT_THROW(m1c({1,3},{1,5}), out_of_range);
    EXPECT_THROW(m1c({4,1},{2,2}), invalid_argument);
    EXPECT_THROW(m1c({1,3},{5,3}), invalid_argument);
    EXPECT_EQ(m1c({1,1},{3,3}), Matrix<complex<double>>(1,1,{17.0+1i}));
    EXPECT_EQ(m1c({1,3},{3,3}), Matrix<complex<double>>(3,1,{17.0+1i,21.0+1i,25.0+1i}));
    EXPECT_EQ(m1c({3,3},{1,2}), Matrix<complex<double>>(1,2,{23.0+1i,24.0+1i}));
    EXPECT_EQ(m1c({1,2},{1,2}), Matrix<complex<double>>(2,2,{15.0+1i,16.0+1i,19.0+1i,20.0+1i}));
    EXPECT_EQ(m1c({1,2},{0,3}), Matrix<complex<double>>(2,4,{14.0+1i,15.0+1i,16.0+1i,17.0+1i,18.0+1i,19.0+1i,20.0+1i,21.0+1i}));
}

TEST(Matrix, setter) {
    Matrix m1(4,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25}
    );

    // vector full set
    Matrix m2 = m1;
    std::vector<double> v = {0,1,2,3,
                             4,5,6,7,
                             8,9,10,11,
                             12,13,14,15};
    m2.set(v);
    for(int i=0; i<4;++i) for(int j=0;j<4;++j) EXPECT_EQ(m2(i,j), v[j+i*4]);

    // vector partial set with double range
    m2 = m1;
    v = {0,1,2,
         3,4,5,
         6,7,8};
    EXPECT_THROW(m2.set({1,4},{1,3}, v), out_of_range);
    EXPECT_THROW(m2.set({1,2},{1,4}, v), out_of_range);
    EXPECT_THROW(m2.set({2,1},{1,3}, v), invalid_argument);
    EXPECT_THROW(m2.set({1,1},{3,1}, v), invalid_argument);
    EXPECT_THROW(m2.set({0,3},{1,3}, v), out_of_range);
    EXPECT_NO_THROW(m2.set({1,3},{1,3}, v));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i==0 || j==0) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), v[j-1 + (i-1)*3]);
    }
    // with only partial values of v
    m2 = m1;
    EXPECT_NO_THROW(m2.set({1,2},{1,2}, v));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i==0 || j==0 || i==3 || j==3) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), v[j-1 + (i-1)*2]);
    }
    m2 = m1;
    EXPECT_NO_THROW(m2.set({0,2},{2,3}, v));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i>2 || j<2) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), v[j-2 + (i)*2]);
    }

    // matrix partial set with double range
    Matrix m3(3,3,v);
    m2 = m1;
    EXPECT_THROW(m2.set({1,4},{1,3}, m3), out_of_range);
    EXPECT_THROW(m2.set({1,2},{1,4}, m3), out_of_range);
    EXPECT_THROW(m2.set({2,1},{1,3}, m3), invalid_argument);
    EXPECT_THROW(m2.set({1,1},{3,1}, m3), invalid_argument);
    EXPECT_THROW(m2.set({0,3},{1,1}, m3), invalid_argument);
    EXPECT_THROW(m2.set(ALL,{1,3}, m3), invalid_argument);
    EXPECT_NO_THROW(m2.set({0,2},{1,3}, m3));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i>2 || j<1) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), m3(i, j-1));
    }

    // matrix partial set with range on columns
    m2 = m1;
    m3 = Matrix(1,3,v);
    EXPECT_THROW(m2.set(4,{1,3}, m3), out_of_range);
    EXPECT_THROW(m2.set(2,{1,4}, m3), out_of_range);
    EXPECT_THROW(m2.set(2,{3,2}, m3), invalid_argument);
    EXPECT_THROW(m2.set(2,{0,3}, m3), invalid_argument);
    EXPECT_NO_THROW(m2.set(2,{1,3}, m3));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i!=2 || j<1) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), m3(0, j-1));
    }

    // matrix partial set with range on rows
    m2 = m1;
    m3 = Matrix(3,1,v);
    EXPECT_THROW(m2.set({1,3}, 4, m3), out_of_range);
    EXPECT_THROW(m2.set({1,4}, 2, m3), out_of_range);
    EXPECT_THROW(m2.set({3,2}, 2, m3), invalid_argument);
    EXPECT_THROW(m2.set({0,3}, 2, m3), invalid_argument);
    EXPECT_NO_THROW(m2.set({1,3}, 2, m3));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(j!=2 || i<1) EXPECT_EQ(m2(i,j), m1(i,j));
        else EXPECT_EQ(m2(i,j), m3(i-1, 0));
    }

    // -- coplex matrix -- //
    Matrix<complex<double>> m1c(4,4,
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
         14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,
         22.0+1i,23.0+1i,24.0+1i,25.0+1i}
    );

    // vector full set
    Matrix m2c(m1c);
    std::vector<complex<double>> v_comp = {
         0.0+1i, 1.0+1i, 2.0+1i, 3.0+1i,
         4.0+1i, 5.0+1i, 6.0+1i, 7.0+1i,
         8.0+1i, 9.0+1i,10.0+1i,11.0+1i,
        12.0+1i,13.0+1i,14.0+1i,15.0+1i
    };
    m2c.set(v_comp);
    for(int i=0; i<4; ++i) for(int j=0;j<4; ++j) EXPECT_EQ(m2c(i,j), v_comp[j+i*4]);

    // vector partial set with double range
    m2c = m1c;
    v_comp = {
        0.0+1i,1.0+1i,2.0+1i,
        3.0+1i,4.0+1i,5.0+1i,
        6.0+1i,7.0+1i,8.0+1i
    };
    EXPECT_THROW(m2c.set({1,4},{1,3}, v_comp), out_of_range);
    EXPECT_THROW(m2c.set({1,2},{1,4}, v_comp), out_of_range);
    EXPECT_THROW(m2c.set({2,1},{1,3}, v_comp), invalid_argument);
    EXPECT_THROW(m2c.set({1,1},{3,1}, v_comp), invalid_argument);
    EXPECT_THROW(m2c.set({0,3},{1,3}, v_comp), out_of_range);
    EXPECT_NO_THROW(m2c.set({1,3},{1,3}, v_comp));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i==0 || j==0) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), v_comp[j-1 + (i-1)*3]);
    }
    // with only partial values of v
    m2c = m1c;
    EXPECT_NO_THROW(m2c.set({1,2},{1,2}, v_comp));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i==0 || j==0 || i==3 || j==3) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), v_comp[j-1 + (i-1)*2]);
    }
    m2c = m1c;
    EXPECT_NO_THROW(m2c.set({0,2},{2,3}, v_comp));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i>2 || j<2) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), v_comp[j-2 + (i)*2]);
    }

    // matrix partial set with double range
    Matrix m3c(3,3,v_comp);
    m2c = m1c;
    EXPECT_THROW(m2c.set({1,4},{1,3}, m3c), out_of_range);
    EXPECT_THROW(m2c.set({1,2},{1,4}, m3c), out_of_range);
    EXPECT_THROW(m2c.set({2,1},{1,3}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set({1,1},{3,1}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set({0,3},{1,1}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set(ALL,{1,3},   m3c), invalid_argument);
    EXPECT_NO_THROW(m2c.set({0,2},{1,3}, m3c));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i>2 || j<1) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), m3c(i, j-1));
    }

    // matrix partial set with range on columns
    m2c = m1c;
    m3c = Matrix(1,3,v_comp);
    EXPECT_THROW(m2c.set(4,{1,3}, m3c), out_of_range);
    EXPECT_THROW(m2c.set(2,{1,4}, m3c), out_of_range);
    EXPECT_THROW(m2c.set(2,{3,2}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set(2,{0,3}, m3c), invalid_argument);
    EXPECT_NO_THROW(m2c.set(2,{1,3}, m3c));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i!=2 || j<1) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), m3c(0, j-1));
    }

    // matrix partial set with range on rows
    m2c = m1c;
    m3c = Matrix(3,1,v_comp);
    EXPECT_THROW(m2c.set({1,3}, 4, m3c), out_of_range);
    EXPECT_THROW(m2c.set({1,4}, 2, m3c), out_of_range);
    EXPECT_THROW(m2c.set({3,2}, 2, m3c), invalid_argument);
    EXPECT_THROW(m2c.set({0,3}, 2, m3c), invalid_argument);
    EXPECT_NO_THROW(m2c.set({1,3}, 2, m3c));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(j!=2 || i<1) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), m3c(i-1, 0));
    }


    // -- from double to complex -- //

    // vector full set
    m2c = m1c;
    v = {0,1,2,3,
         4,5,6,7,
         8,9,10,11,
         12,13,14,15};

    m2c.set(v);
    for(int i=0; i<4; ++i) for(int j=0;j<4; ++j) EXPECT_EQ(m2c(i,j), v[j+i*4]);

    // vector partial set with double range
    m2c = m1c;
    v = {0,1,2,
         3,4,5,
         6,7,8};
    EXPECT_THROW(m2c.set({1,4},{1,3}, v), out_of_range);
    EXPECT_THROW(m2c.set({1,2},{1,4}, v), out_of_range);
    EXPECT_THROW(m2c.set({2,1},{1,3}, v), invalid_argument);
    EXPECT_THROW(m2c.set({1,1},{3,1}, v), invalid_argument);
    EXPECT_THROW(m2c.set({0,3},{1,3}, v), out_of_range);
    EXPECT_NO_THROW(m2c.set({1,3},{1,3}, v));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i==0 || j==0) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), v[j-1 + (i-1)*3]);
    }
    // with only partial values of v
    m2c = m1c;
    EXPECT_NO_THROW(m2c.set({1,2},{1,2}, v));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i==0 || j==0 || i==3 || j==3) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), v[j-1 + (i-1)*2]);
    }
    m2c = m1c;
    EXPECT_NO_THROW(m2c.set({0,2},{2,3}, v));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i>2 || j<2) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), v[j-2 + (i)*2]);
    }

    // matrix partial set with double range
    m2c = m1c;
    m3c = Matrix(3,3,v);
    EXPECT_THROW(m2c.set({1,4},{1,3}, m3c), out_of_range);
    EXPECT_THROW(m2c.set({1,2},{1,4}, m3c), out_of_range);
    EXPECT_THROW(m2c.set({2,1},{1,3}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set({1,1},{3,1}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set({0,3},{1,1}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set(ALL,{1,3},   m3c), invalid_argument);
    EXPECT_NO_THROW(m2c.set({0,2},{1,3}, m3c));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i>2 || j<1) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), m3c(i, j-1));
    }

    // matrix partial set with range on columns
    m2c = m1c;
    m3c = Matrix(1,3,v);
    EXPECT_THROW(m2c.set(4,{1,3}, m3c), out_of_range);
    EXPECT_THROW(m2c.set(2,{1,4}, m3c), out_of_range);
    EXPECT_THROW(m2c.set(2,{3,2}, m3c), invalid_argument);
    EXPECT_THROW(m2c.set(2,{0,3}, m3c), invalid_argument);
    EXPECT_NO_THROW(m2c.set(2,{1,3}, m3c));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(i!=2 || j<1) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), m3c(0, j-1));
    }

    // matrix partial set with range on rows
    m2c = m1c;
    m3c = Matrix(3,1,v);
    EXPECT_THROW(m2c.set({1,3}, 4, m3c), out_of_range);
    EXPECT_THROW(m2c.set({1,4}, 2, m3c), out_of_range);
    EXPECT_THROW(m2c.set({3,2}, 2, m3c), invalid_argument);
    EXPECT_THROW(m2c.set({0,3}, 2, m3c), invalid_argument);
    EXPECT_NO_THROW(m2c.set({1,3}, 2, m3c));
    for(int i=0; i<4; ++i) for(int j=0;j<4;++j){
        if(j!=2 || i<1) EXPECT_EQ(m2c(i,j), m1c(i,j));
        else EXPECT_EQ(m2c(i,j), m3c(i-1, 0));
    }
}

TEST(Matrix, reshape){
    Matrix m1(4,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25}
    );
    Matrix m2(2,8,
        {10,11,12,13,14,15,16,17,
         18,19,20,21,22,23,24,25}
    );

    EXPECT_THROW(Matrix().reshape(0,0), runtime_error);
    EXPECT_THROW(m1.reshape(3,5), invalid_argument);
    EXPECT_EQ(m1.reshape(2,8), m2);

    // -- complex matrix -- //
    Matrix<complex<double>> m1c(4,4,
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
         14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,
         22.0+1i,23.0+1i,24.0+1i,25.0+1i}
    );
    Matrix<complex<double>> m2c(2,8,
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,22.0+1i,23.0+1i,24.0+1i,25.0+1i}
    );

    EXPECT_THROW(Matrix().reshape(0,0), runtime_error);
    EXPECT_THROW(m1c.reshape(3,5), invalid_argument);
    EXPECT_EQ(m1c.reshape(2,8), m2c);
    
}

TEST(Matrix, swap){
    Matrix m1 = RandMat(3,4);
    Matrix m2 = RandMat(5,6);

    // make a copy
    Matrix copy1(m1);
    Matrix copy2(m2);

    //swap
    EXPECT_NO_THROW(swap(m1, m2));
    EXPECT_EQ(copy1, m2);
    EXPECT_EQ(copy2, m1);


    // -- complex matrix -- //
    Matrix m1c = RandMat<complex<double>>(3,4);
    Matrix m2c = RandMat<complex<double>>(5,6);

    // make a copy
    Matrix copy1c(m1c);
    Matrix copy2c(m2c);

    //swap
    EXPECT_NO_THROW(swap(m1c, m2c));
    EXPECT_EQ(copy1c, m2c);
    EXPECT_EQ(copy2c, m1c);

    // swap between different matrices
    // compilation error: not defined
    // swap(m1c, m2);
}

TEST(Matrix, swap_rows){
    Matrix m1(4,4,
        {10,11,12,13,
         14,15,16,17,
         18,19,20,21,
         22,23,24,25}
    );

    // throws out of range
    EXPECT_THROW(m1.swap_rows(1,4), out_of_range);
    EXPECT_THROW(m1.swap_rows(4,1), out_of_range);

    // swap
    m1.swap_rows(1,3);
    EXPECT_EQ(m1, Matrix(4,4,
        {10,11,12,13,
         22,23,24,25,
         18,19,20,21,
         14,15,16,17}
    ));
}

TEST(Matrix, swap_columns){
    Matrix<complex<double>> m1(4,4,
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
         14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,
         22.0+1i,23.0+1i,24.0+1i,25.0+1i}
    );

    // throws out of range
    EXPECT_THROW(m1.swap_cols(1,4), out_of_range);
    EXPECT_THROW(m1.swap_cols(4,1), out_of_range);

    m1.swap_cols(0,2);
    EXPECT_EQ(m1, Matrix<complex<double>>(4,4,
        {12.0+1i,11.0+1i,10.0+1i,13.0+1i,
         16.0+1i,15.0+1i,14.0+1i,17.0+1i,
         20.0+1i,19.0+1i,18.0+1i,21.0+1i,
         24.0+1i,23.0+1i,22.0+1i,25.0+1i}
    ));
}

TEST(Matrix, transpose){
//     Matrix m2;

//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_NO_THROW(m2 = m1.t());
//     EXPECT_EQ(m2.r(), (uint)4);
//     EXPECT_EQ(m2.c(), (uint)3);
//     for(int i=0; i<3; ++i) for(int j=0;j<4;++j){
//         EXPECT_EQ(m1(i,j), m2(j,i));
//     }

//     m1 = Matrix(1,5,
//         {10,11,12,13,14}
//     );
//     EXPECT_NO_THROW(m2 = m1.t());
//     EXPECT_EQ(m2.r(), (uint)5);
//     EXPECT_EQ(m2.c(), (uint)1);
//     for(int i=0; i<1; ++i) for(int j=0;j<5;++j){
//         EXPECT_EQ(m1(i,j), m2(j,i));
//     }
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
    EXPECT_TRUE( Matrix(2, 4, v1) == Matrix(2, 4, v1));
    EXPECT_FALSE(Matrix(3, 4, v1) == Matrix(3, 4, v2));

    // -- complex matrix -- //
    std::vector<complex<double>> v1c = 
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
         14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,
         22.0+1i,23.0+1i,24.0+1i,25.0+1i};
    std::vector<complex<double>> v2c = 
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
          4.0+1i, 5.0+1i, 6.0+1i, 7.0+1i,
          8.0+1i, 9.0+1i,10.0+1i,11.0+1i,
         12.0+1i,13.0+1i,14.0+1i,15.0+1i};

    EXPECT_FALSE(Matrix(3, 4, v1c) == Matrix(4, 3, v1c));
    EXPECT_FALSE(Matrix(3, 4, v1c) == Matrix(3, 3, v1c));
    EXPECT_TRUE( Matrix(2, 4, v1c) == Matrix(2, 4, v1c));
    EXPECT_FALSE(Matrix(3, 4, v1c) == Matrix(3, 4, v2c));

    // -- double and complex
    EXPECT_FALSE(Matrix<complex<double>>(3, 4, v1) == Matrix<double>(4, 3, v1));
    EXPECT_FALSE(Matrix<complex<double>>(3, 4, v1) == Matrix<double>(3, 3, v1));
    EXPECT_TRUE( Matrix<complex<double>>(2, 4, v1) == Matrix<double>(2, 4, v1));
    EXPECT_FALSE(Matrix<complex<double>>(3, 4, v1) == Matrix<double>(3, 4, v2));
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
    EXPECT_TRUE(Matrix(3, 4, v1) != Matrix(3, 4, v2));

    // -- complex matrix -- //
    std::vector<complex<double>> v1c = 
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
         14.0+1i,15.0+1i,16.0+1i,17.0+1i,
         18.0+1i,19.0+1i,20.0+1i,21.0+1i,
         22.0+1i,23.0+1i,24.0+1i,25.0+1i};
    std::vector<complex<double>> v2c = 
        {10.0+1i,11.0+1i,12.0+1i,13.0+1i,
          4.0+1i, 5.0+1i, 6.0+1i, 7.0+1i,
          8.0+1i, 9.0+1i,10.0+1i,11.0+1i,
         12.0+1i,13.0+1i,14.0+1i,15.0+1i};

    EXPECT_TRUE(Matrix(3, 4, v1c) != Matrix(4, 3, v1c));
    EXPECT_TRUE(Matrix(3, 4, v1c) != Matrix(3, 3, v1c));
    EXPECT_FALSE( Matrix(2, 4, v1c) != Matrix(2, 4, v1c));
    EXPECT_TRUE(Matrix(3, 4, v1c) != Matrix(3, 4, v2c));

    // -- double and complex
    EXPECT_TRUE(Matrix<complex<double>>(3, 4, v1) != Matrix<double>(4, 3, v1));
    EXPECT_TRUE(Matrix<complex<double>>(3, 4, v1) != Matrix<double>(3, 3, v1));
    EXPECT_FALSE( Matrix<complex<double>>(2, 4, v1) != Matrix<double>(2, 4, v1));
    EXPECT_TRUE(Matrix<complex<double>>(3, 4, v1) != Matrix<double>(3, 4, v2));
}

TEST(Matrix, sum_operator){
//     std::vector<double> v1 =  {10,11,12,13,14,15,16,17,18,19,20,21};
//     std::vector<double> v2 = v1;
//     Matrix m1(3,4,v1);

//     // int
//     for(int i=0; i<v2.size(); ++i) v2[i]+=7;
//     EXPECT_EQ(m1+7, Matrix(3,4, v2));
//     m1+=7;
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)+=3);
//     EXPECT_NO_THROW(Matrix(3,0)+3);

//     // double
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     for(int i=0; i<v2.size(); ++i) v2[i]+=(-3.14);
//     EXPECT_EQ(m1+(-3.14) ,Matrix(3,4, v2));
//     m1+=(-3.14);
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)+=7.72);
//     EXPECT_NO_THROW(Matrix(3,0)+8.13);

//     // matrices
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     Matrix m2 = m1 + 2.28;
//     for(int i=0; i<v2.size(); ++i) v2[i] = v2[i] + v2[i] + 2.28;
//     EXPECT_EQ((m1+m2), Matrix(3,4, v2));
//     m1+=m2;
//     EXPECT_EQ(m1, Matrix(3,4, v2));

//     Matrix m3(4,3, v2);
//     EXPECT_THROW(m1+m3, invalid_argument);
//     EXPECT_THROW(m1+=m3, invalid_argument);
}

TEST(Matrix, subtract_operator){
//     std::vector<double> v1 =  {10,11,12,13,14,15,16,17,18,19,20,21};
//     std::vector<double> v2 = v1;
//     Matrix m1(3,4,v1);

//     // int
//     for(int i=0; i<v2.size(); ++i) v2[i]-=7;
//     EXPECT_EQ(m1-7, Matrix(3,4, v2));
//     m1-=7;
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)-=3);
//     EXPECT_NO_THROW(Matrix(3,0)-3);

//     // double
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     for(int i=0; i<v2.size(); ++i) v2[i]-=(-3.14);
//     EXPECT_EQ(m1-(-3.14), Matrix(3,4, v2));
//     m1-=(-3.14);
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)-=7.72);
//     EXPECT_NO_THROW(Matrix(3,0)-8.13);

//     // matrices
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     Matrix m2 = m1 + 2;
//     for(int i=0; i<v2.size(); ++i) v2[i] = v2[i] - v2[i] - 2;
//     EXPECT_EQ((m1-m2), Matrix(3,4, v2));
//     // cout << m1-m2 << endl;
//     // cout << Matrix(3,4, v2) << endl;
//     m1-=m2;
//     EXPECT_EQ(m1, Matrix(3,4, v2));

//     Matrix m3(4,3, v2);
//     EXPECT_THROW(m1-m3, invalid_argument);
//     EXPECT_THROW(m1-=m3, invalid_argument);

}

TEST(Matrix, multiply_operator){
//     std::vector<double> v1 =  {1,3,5,9,1,3,1,7,4,3,9,7};
//     std::vector<double> v2 = v1;
//     Matrix m1(3,4,v1);

//     // int
//     for(int i=0; i<v2.size(); ++i) v2[i]*=2;
//     EXPECT_EQ(m1*2, Matrix(3,4, v2));
//     m1*=2;
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)*=2);
//     EXPECT_NO_THROW(Matrix(3,0)*2);

//     // double
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     for(int i=0; i<v2.size(); ++i) v2[i]*=2.5;
//     EXPECT_EQ(m1*2.5, Matrix(3,4, v2));
//     m1*=2.5;
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)*=2.5);
//     EXPECT_NO_THROW(Matrix(3,0)*2.5);

//     // matrices
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     Matrix m2 = m1({1,2}, {0,3}).t();
//     EXPECT_EQ((m1*m2), Matrix(3,2, {78,121,60,71,71,155}));
//     m1*=m2;
//     EXPECT_EQ(m1, Matrix(3,2, {78,121,60,71,71,155}));
//     m1 = Matrix(3,4,v1);
//     EXPECT_EQ(m1(1,{1,3})*m1({0,2},2), Matrix(1,1, {79}));
//     m1*=m1.t();
//     EXPECT_EQ(m1, Matrix(3,3, {116,78,121,78,60,71,121,71,155}));
//     m1 = Matrix(3,4,v1);
//     m2 = m1({1,2}, {0,3});
//     EXPECT_THROW(m1*m2, invalid_argument);
//     EXPECT_THROW(m1*=m1, invalid_argument);
}

TEST(Matrix, concatenate_operators){
//     Matrix m1(4,6,{ 14,   4, -16, -20, -3,  -4,
//                     17,  17, -16,   5,  6, -13,
//                     -2,  12,  -7,  15, -1,  -3,
//                     19,   4,  -3,   5, 16,  -9}); 
    
//     // concatenate per rows
//     EXPECT_EQ(m1,
//         Matrix(1,6, {14,   4, -16, -20, -3,  -4}) |
//         Matrix(1,6, {17,  17, -16,   5,  6, -13}) |
//         Matrix(1,6, {-2,  12,  -7,  15, -1,  -3}) |
//         Matrix(1,6, {19,   4,  -3,   5, 16,  -9})
//     );

//     // concatenate per columns
//     EXPECT_EQ(m1,
//         Matrix(4,1, {14, 17, -2, 19}) &
//         Matrix(4,1, {4, 17, 12, 4}) &
//         Matrix(4,1, {-16, -16, -7, -3}) &
//         Matrix(4,1, {-20, 5, 15, 5}) &
//         Matrix(4,1, {-3, 6, -1, 16}) &
//         Matrix(4,1, { -4, -13, -3, -9})
//     );
}



TEST(Matrix, to_vec){
    Matrix<complex<double>> m1(2,3,{10.0+1i,11.0+1i,13.0+1i,14.0+1i,15.0+1i,17.0+1i});
    Matrix<complex<double>> m2(1,6,{10.0+1i,11.0+1i,13.0+1i,14.0+1i,15.0+1i,17.0+1i});
    Matrix<complex<double>> m3(6,1,{10.0+1i,11.0+1i,13.0+1i,14.0+1i,15.0+1i,17.0+1i});

    EXPECT_THROW(Matrix().to_c_vec(), runtime_error);
    EXPECT_THROW(Matrix().to_r_vec(), runtime_error);
    EXPECT_EQ(m1.to_c_vec(), m3);
    EXPECT_EQ(m1.to_r_vec(), m2);
}

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
    EXPECT_EQ(m1({3,3},1).norm2(), 2);

    // -- complex matrix -- //
    Matrix<complex<double>> m2(4,4,
        {1,     3,5,9,
         1,3.0+4i,1,7,
         4,3.0+4i,9,7,
         5,3.0+4i,0,9}
    );

    EXPECT_THROW(m2.norm2(), invalid_argument);
    EXPECT_EQ   (m2(2,{0,3}).norm2(), sqrt(171));
    EXPECT_EQ   (m2({0,2},1).norm2(), sqrt(59));
    EXPECT_EQ   (m2({3,3},1).norm2(), 5);
}

TEST(Matrix, normalize){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_THROW(m1.normalize(), invalid_argument);

    m1 = Matrix(1,4, {4,4,4,4});
    EXPECT_EQ(m1.normalize(), Matrix(1,4, {0.5,0.5,0.5,0.5}));

    // -- complex matrix -- //
    
}

TEST(Matrix, normalize_self){
//     Matrix m1(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );
//     EXPECT_THROW(m1.normalize_self(), invalid_argument);

//     Matrix m2(1,4, {4,4,4,4});
//     m2.normalize_self();
//     EXPECT_EQ(m2, Matrix(1,4, {0.5,0.5,0.5,0.5}));
}



TEST(Matrix, IdMat){
//     Matrix m1;
    
//     EXPECT_NO_THROW(m1 = IdMat(4));
//     EXPECT_EQ(m1.r(), (uint)4);
//     EXPECT_EQ(m1.c(), (uint)4);

//     for(uint i=0; i<4; ++i){
//         for(uint j=0; j<4; ++j){
//             if(i == j) EXPECT_EQ(m1(i,j), 1);
//             else EXPECT_EQ(m1(i,j), 0);
//         }
//     }
}

TEST(Matrix, Ones){
//     Matrix m1;
    
//     EXPECT_NO_THROW(m1 = Ones(4,3));
//     EXPECT_EQ(m1.r(), (uint)4);
//     EXPECT_EQ(m1.c(), (uint)3);

//     for(uint i=0; i<4; ++i) for(uint j=0; j<3; ++j) EXPECT_EQ(m1(i,j), 1);
}

TEST(Matrix, RandMat){
//     Matrix m1;
    
//     EXPECT_NO_THROW(m1 = RandMat(5,6));
//     EXPECT_EQ(m1.r(), (uint)5);
//     EXPECT_EQ(m1.c(), (uint)6);

//     const double * v = m1.v();
//     double tmp = v[0];
//     for(uint i=1; i<30; ++i){
//         EXPECT_TRUE(tmp != v[i]);
//         tmp = v[i];
//     }
}

TEST(Matrix, diag){
//     Matrix m1;
    
//     EXPECT_NO_THROW(m1 = diag({1,2,3,4}));
//     EXPECT_EQ(m1.r(), (uint)4);
//     EXPECT_EQ(m1.c(), (uint)4);
//     for(uint i=0; i<4; ++i){
//         for(uint j=0; j<4; ++j){
//             if(i == j) EXPECT_EQ(m1(i,j), i + 1);
//             else EXPECT_EQ(m1(i,j), 0);
//         }
//     }

//     Matrix m2(1,4,{1,3,5,9});

//     EXPECT_NO_THROW(m1 = diag(m2));
//     EXPECT_EQ(m1.r(), (uint)4);
//     EXPECT_EQ(m1.c(), (uint)4);
//     for(uint i=0; i<4; ++i){
//         for(uint j=0; j<4; ++j){
//             if(i == j) EXPECT_EQ(m1(i,j), m2(i));
//             else EXPECT_EQ(m1(i,j), 0);
//         }
//     }
}



TEST(Matrix, submat_del){
//     Matrix m2;

//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.submat_del(3,2), invalid_argument);
//     EXPECT_THROW(m1.submat_del(0,4), invalid_argument);
//     EXPECT_NO_THROW(m2 = m1.submat_del(1,2));
//     EXPECT_EQ(m2.r(), (uint)2);
//     EXPECT_EQ(m2.c(), (uint)3);
//     for(int i=0; i<2; ++i) for(int j=0;j<3;++j){
//         uint k = i > 0 ? i+1 : i;
//         uint l = j > 1 ? j+1 : j;
//         EXPECT_EQ(m2(i,j), m1(k,l));
//     }

//     m1 = Matrix(1,5,
//         {10,11,12,13,14}
//     );
//     EXPECT_NO_THROW(m2 = m1.submat_del(0,3));
//     EXPECT_EQ(m2.r(), (uint)0);
//     EXPECT_EQ(m2.c(), (uint)4);
//     EXPECT_THROW(m2(0,0), out_of_range);

}

TEST(Matrix, determinant){
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.det(), invalid_argument);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );
    
//     // 1x1
//     EXPECT_DOUBLE_EQ(m1({1,1},0).det(), 1);
//     EXPECT_DOUBLE_EQ(m1({2,2},{1,1}).det(), 3);
//     // 2x2
//     EXPECT_DOUBLE_EQ(m1({1,2},{2,3}).det(), -56);
//     // 3x3
//     EXPECT_DOUBLE_EQ(m1({1,3},{0,2}).det(), 110);
//     // 4x4
//     EXPECT_DOUBLE_EQ(m1.det(), -376);

//     m1 = Matrix(4,4,
//         {1,3,5,0,
//          1,3,1,0,
//          4,3,9,0,
//          5,2,0,0}
//     );
//     EXPECT_DOUBLE_EQ(m1.det(), 0);

}

TEST(Matrix, minor){
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.minor(1,1), invalid_argument);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );
//     double minor;
//     EXPECT_NO_THROW(minor = m1.minor(1,1));
//     EXPECT_DOUBLE_EQ(minor, -329.0);
//     EXPECT_NO_THROW(minor = m1.minor(2,1));
//     EXPECT_DOUBLE_EQ(minor, 94.0);
}

TEST(Matrix, cofactor){
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.cof(1,1), invalid_argument);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );
//     double cofactor;
//     EXPECT_NO_THROW(cofactor = m1.cof(1,1));
//     EXPECT_DOUBLE_EQ(cofactor, -329.0);
//     EXPECT_NO_THROW(cofactor = m1.cof(2,1));
//     EXPECT_DOUBLE_EQ(cofactor, -94.0);
}

TEST(Matrix, cofactor_matrix){
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.cof_mat(), invalid_argument);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );

//     Matrix m2;
//     EXPECT_NO_THROW(m2 = m1.cof_mat());
//     EXPECT_EQ(m2, Matrix(4,4,{
//         104, 235, -39, -110,
//         -16, -329, 53,   82,
//         -56, -94, -26,   52,
//         -48, 94,   18,  -36}));
    
}

TEST(Matrix, adj){
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.adj(), invalid_argument);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );

//     Matrix m2;
//     EXPECT_NO_THROW(m2 = m1.adj());
//     EXPECT_EQ(m2, Matrix(4,4,{
//         104, 235, -39, -110,
//         -16, -329, 53,   82,
//         -56, -94, -26,   52,
//         -48, 94,   18,  -36}).t());
}

TEST(Matrix, inv){
//     // not square
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.inv(), invalid_argument);

//     // singular
//     m1 = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,0,0,0}
//     );
//     EXPECT_THROW(m1.inv(), runtime_error);
//     m1 = Matrix (4,4,
//         {1,3, 4, 0,
//          0,2, 1, 0,
//          2,1,-3, 0,
//          7,0, 2, 0}
//     );
//     EXPECT_THROW(m1.inv(), runtime_error);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );
//     Matrix res(4,4, 
//         {-104,  16,  56,  48,
//          -235, 329,  94, -94,
//            39, -53,  26, -18,
//           110, -82, -52,  36}
//     );
//     EXPECT_EQ(m1.inv(), res/376);
//     EXPECT_EQ(m1*m1.inv(), IdMat(4));
}

TEST(Matrix, pinv){
//     Matrix m1(4,3,
//         {1,3,4,
//          2,0,2,
//          1,-2,2,
//          1,-3,2}
//     );
//     Matrix m2;

//     EXPECT_NO_THROW(m2 = m1.pinv_left());
//     EXPECT_EQ(m2*m1, IdMat(3));

//     m1 = Matrix(3,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2});

//     EXPECT_NO_THROW(m2 = m1.pinv_right());
//     EXPECT_EQ(m1*m2, IdMat(3));
}



TEST(Matrix, is_sing){
//     Matrix m1(3,4,
//         {10,11,12,13,
//          14,15,16,17,
//          18,19,20,21}
//     );
//     EXPECT_THROW(m1.det(), invalid_argument);

//     m1 = Matrix(4,4,
//         {1,3,5,9,
//          1,3,1,7,
//          4,3,9,7,
//          5,2,0,9}
//     );
//     EXPECT_FALSE(m1.is_sing());

//     m1 = Matrix(3,3,
//         {1,1,0,
//          1,2,1,
//          1,3,2}
//     );
//     EXPECT_TRUE(m1.is_sing());
}

TEST(Matrix, is_vec){
    EXPECT_FALSE(
        Matrix<complex<double>>(2,3,
            {10.0+1i,11,13,
            14,15,17.0+1i}
        ).is_vec()
    );
    EXPECT_TRUE(
        Matrix(1,3,
            {10,11,13}
        ).is_vec()
    );
    EXPECT_TRUE(
        Matrix<complex<double>>(3,1,
            {10,11.0+1i,13}
        ).is_vec()
    );
}

TEST(Matrix, is_upper_triang){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_FALSE(m1.is_upper_triang());

    Matrix<complex<double>> m2(4,4,
        {1,3,5,9,
         0,3,1,7,
         0.0+0i,0,9,7,
         0,0.0+0i,0.0+0i,9}
    );
    EXPECT_TRUE(m2.is_upper_triang());

    m2 = Matrix<complex<double>>(4,4,
        {1,3,5,9,
         0,3,1,7,
         0.0+1i,0,9,7,
         0,0.0+1i,0.0+1i,9}
    );
    EXPECT_FALSE(m2.is_upper_triang());
}

TEST(Matrix, is_lower_triang){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_FALSE(m1.is_lower_triang());

    Matrix<complex<double>> m2(4,4,
        {1,0,0,0.0+0i,
         1,3,0.0+0i,0,
         4,3,9,0.0+0i,
         5,2,0,9}
    );
    EXPECT_TRUE(m2.is_lower_triang());

    m2 = Matrix<complex<double>>(4,4,
        {1,0,0,0,
         1,3,0,0,
         4,3,9,0.0+1i,
         5,2,0,9}
    );
    EXPECT_FALSE(m2.is_lower_triang());
}

TEST(Matrix, is_upper_hessenberg){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         0,0,0,9}
    );
    EXPECT_FALSE(m1.is_upper_hessenberg());

    Matrix<complex<double>> m2(4,4,
        {1,3,5,9,
         1,3.0+1i,1,7,
         0,1,9,7.0+1i,
         0,0.0+0i,1,9}
    );
    EXPECT_TRUE(m2.is_upper_hessenberg());

    m2 = Matrix<complex<double>>(4,4,
        {1.0+1i,3,5,9,
         1,3,1.0+1i,7,
         0,1,9,7,
         0,0.0+0.1i,1,9}
    );
    EXPECT_FALSE(m2.is_upper_hessenberg());
}

TEST(Matrix, is_lower_hessenberg){
    Matrix m1(4,4,
        {1,3,5,9,
         1,3,1,7,
         4,3,9,7,
         5,2,0,9}
    );
    EXPECT_FALSE(m1.is_lower_hessenberg());

    m1 = Matrix(4,4,
        {1,1,0,0,
         1,3,1,0,
         4,3,9,1,
         5,2,0,9}
    );
    EXPECT_TRUE(m1.is_lower_hessenberg());

    Matrix<complex<double>> m2(4,4,
        {1,1,0,0,
         1,3,1,0.0+0.1i,
         4.0+1i,3,9,1,
         5,2,0.0+1i,9}
    );
    EXPECT_FALSE(m2.is_lower_hessenberg());
}



// TEST(Matrix, qr_dec){
//     Matrix Q, R;

//     Matrix A(3,3, {2, -2, 18,
//                    2, 1, 0,
//                    1, 2, 0});
//     EXPECT_NO_THROW(A.qr_dec(Q,R));
//     EXPECT_EQ(Q.t()*Q, IdMat(3));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A, Q*R);

//     A = Matrix (3,3, {12, -51, 4,
//                       6, 167, -68,
//                       -4, 24, -41});
//     EXPECT_NO_THROW(A.qr_dec(Q,R));
//     EXPECT_EQ(Q.t()*Q, IdMat(3));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A, Q*R);

//     A = Matrix(4,6,{ 14,   4, -16, -20, -3,  -4,
//                      17,  17, -16,   5,  6, -13,
//                      -2,  12,  -7,  15, -1,  -3,
//                      19,   4,  -3,   5, 16,  -9});
//     EXPECT_NO_THROW(A.qr_dec(Q,R));
//     EXPECT_EQ(Q.t()*Q, IdMat(4));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A, Q*R);
// }

// TEST(Matrix, qrp_dec){
//     Matrix Q, R, P;

//     Matrix A(3,3, {2, -2, 18,
//                    2, 1, 0,
//                    1, 2, 0});
//     EXPECT_NO_THROW(A.qrp_dec(Q,R,P));
//     EXPECT_EQ(Q.t()*Q, IdMat(3));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A*P, Q*R);

//     A = Matrix (3,3, {12, -51, 4,
//                       6, 167, -68,
//                       -4, 24, -41});
//     EXPECT_NO_THROW(A.qrp_dec(Q,R,P));
//     EXPECT_EQ(Q.t()*Q, IdMat(3));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A*P, Q*R);

//     A = Matrix(4,6,{ 14,   4, -16, -20, -3,  -4,
//                      17,  17, -16,   5,  6, -13,
//                      -2,  12,  -7,  15, -1,  -3,
//                      19,   4,  -3,   5, 16,  -9});
//     EXPECT_NO_THROW(A.qrp_dec(Q,R,P));
//     EXPECT_EQ(Q.t()*Q, IdMat(4));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A*P, Q*R);

//     A = Matrix(6,4,{ 14,   4, -16, -20, -3,  -4,
//                      17,  17, -16,   5,  6, -13,
//                      -2,  12,  -7,  15, -1,  -3,
//                      19,   4,  -3,   5, 16,  -9});
//     EXPECT_NO_THROW(A.qrp_dec(Q,R,P));
//     EXPECT_EQ(Q.t()*Q, IdMat(6));
//     EXPECT_TRUE(R.is_upper_triang());
//     EXPECT_EQ(A*P, Q*R);
// }

// TEST(Matrix, lup_dec){
//     // matrix not square
//     Matrix m1,L,U,P;
//     m1 = Matrix (3,4,
//         {1,3,5,9,
//          0,2,1,7,
//          4,1,8,2,
//          5,2,1,9}
//     );
//     EXPECT_THROW(m1.lup_dec(L,U,P), invalid_argument);

//     // matrix decomposable
//     m1 = Matrix (4,4,
//         {1,3,5,2,
//          0,2,1,6,
//          4,1,3,2,
//          5,2,1,4}
//     );
//     EXPECT_NO_THROW(m1.lup_dec(L,U,P));
//     EXPECT_EQ(L*U, P*m1);

//     m1 = Matrix (4,4,
//         {2, 0, 2, 0.6,
//          3, 3, 4, -2,
//          5, 5, 4, 2,
//         -1, -2, 3.4, -1}
//     );
//     EXPECT_NO_THROW(m1.lup_dec(L,U,P));
//     EXPECT_EQ(L*U, P*m1);

// }

// TEST(Matrix, hessenberg_dec){
//     Matrix m1, Q, H;

//     // matrix not square
//     m1 = Matrix (3,4,
//         {1,3,5,9,
//          0,2,1,7,
//          4,1,8,2,
//          5,2,1,9}
//     );
//     EXPECT_THROW(m1.hessenberg_dec(Q, H), invalid_argument);

//     // matrix decomposable
//     m1 = Matrix (4,4,
//         {7.5231e-01,   8.7419e-01,   3.6122e-01,   4.6593e-01,
//          6.4349e-01,   2.9453e-01,   4.3203e-01,   8.8371e-03,
//          1.4175e-01,   8.3325e-01,   6.4892e-01,   8.3927e-02,
//          8.1433e-01,   9.5796e-01,   9.0255e-01,   1.0307e-01}
//     );
//     EXPECT_NO_THROW(m1.hessenberg_dec(Q, H));
//     EXPECT_EQ(Q * Q.t(), IdMat(4));
//     EXPECT_TRUE(H.is_upper_hessenberg());
//     EXPECT_EQ(Q * H * Q.t(), m1);
    
// }



// TEST(Matrix, backward_sub){
//     Matrix U,b;
//     U = Matrix(4,4,
//         {1,3,5,9,
//          0,1,1,7,
//          0,0,2,7,
//          0,0,0,3}
//     );
//     b = Matrix(4,1, {2,6,1,3});

//     // U not square
//     EXPECT_THROW(Matrix<T>::backward_sub(U({0,3},{1,3}), b), invalid_argument);
//     // b.rows != U.cols
//     EXPECT_THROW(Matrix<T>::backward_sub(U, b({1,3},0)), invalid_argument);
//     // underdetermined
//     U(2,2) = 0;
//     EXPECT_THROW(Matrix<T>::backward_sub(U, b), runtime_error);
//     U(2,2) = 2;
    
//     // solve
//     // cout << Matrix<T>::backward_sub(U, b) << endl;
//     EXPECT_EQ(Matrix<T>::backward_sub(U, b), 
//         Matrix(4,1, {2,2,-3,1})
//     );

// }

// TEST(Matrix, forward_sub){
//     Matrix L,b;
//     L = Matrix(4,4,
//         {1,3,5,9,
//          0,1,1,7,
//          0,0,2,7,
//          0,0,0,3}
//     ).t();
//     b = Matrix(4,1, {2,6,1,3});

//     // L not square
//     EXPECT_THROW(Matrix<T>::forward_sub(L({0,3},{1,3}), b), invalid_argument);
//     // b.rows != L.cols
//     EXPECT_THROW(Matrix<T>::forward_sub(L, b({1,3},0)), invalid_argument);
//     // underdetermined
//     L(2,2) = 0;
//     EXPECT_THROW(Matrix<T>::forward_sub(L, b), runtime_error);
//     L(2,2) = 2;
    
//     // solve
//     // cout << Matrix<T>::forward_sub(L, b) << endl;
//     EXPECT_EQ(Matrix<T>::forward_sub(L, b), 
//         Matrix(4,1, {2,0,-4.5,5.5})
//     );

// }

// TEST(Matrix, matrix_l_divide){
//     Matrix A,b,C;
//     A = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,2,1,-1}
//     );
//     b = Matrix(4,1, {2,4,1,3});

//     // shape doesn't match
//     EXPECT_THROW(Matrix<T>::matrix_l_divide(A, b({1,3},0)), invalid_argument);
    
//     // underdetermined
//     C = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,0,0,0}
//     );
//     EXPECT_THROW(Matrix<T>::matrix_l_divide(C, b), runtime_error);
//     C = Matrix (4,4,
//         {1,3, 4, 0,
//          0,2, 1, 0,
//          2,1,-3, 0,
//          7,0, 2, 0}
//     );
//     EXPECT_THROW(Matrix<T>::matrix_l_divide(C, b), runtime_error);

//     // solution
//     EXPECT_NO_THROW(C = Matrix<T>::matrix_l_divide(A,b));
//     EXPECT_EQ(A*C, b);

//     // overdetermined
//     Matrix<T>::set_double_precision(8);
//     EXPECT_NO_THROW(C = Matrix<T>::matrix_l_divide(A(ALL,{1,3}), b));
//     EXPECT_EQ(C, 
//         Matrix(3,1,{1.42860766, 
//                     -0.26654831, 
//                     -0.62490489})
//     );
//     Matrix<T>::set_double_precision();
// }

// TEST(Matrix, divide_operator){
//     std::vector<double> v1 =  {1,3,5,9,1,3,1,7,4,3,9,7};
//     std::vector<double> v2 = v1;
//     Matrix m1(3,4,v1);

//     // int
//     for(int i=0; i<v2.size(); ++i) v2[i]/=2;
//     EXPECT_EQ(m1/2, Matrix(3,4, v2));
//     m1/=2;
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)/=2);
//     EXPECT_NO_THROW(Matrix(3,0)/2);

//     // double
//     m1 = Matrix(3,4,v1);
//     v2 = v1;
//     for(int i=0; i<v2.size(); ++i) v2[i]/=2.5;
//     EXPECT_EQ(m1/2.5, Matrix(3,4, v2));
//     m1/=2.5;
//     EXPECT_EQ(m1, Matrix(3,4, v2));
//     EXPECT_NO_THROW(Matrix(0,3)/=2.5);
//     EXPECT_NO_THROW(Matrix(3,0)/2.5);

//     // matrices
//     Matrix A(4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,2,1,-1}
//     );
//     Matrix b(1,4, {2,4,1,3});
//     // shape doesn't match
//     EXPECT_THROW((b(0,{1,3}) / A), invalid_argument);
//     // underdetermined
//     Matrix C(4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,0,0,0}
//     );
//     EXPECT_THROW(b / C, runtime_error);
//     // solution
//     EXPECT_NO_THROW(C = b/A);
//     EXPECT_EQ(C*A, b);
//     // shape doesn't match
//     EXPECT_THROW((b(0,{1,3})/=A), invalid_argument);
//     // underdetermined
//     C = Matrix(4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,0,0,0}
//     );
//     EXPECT_THROW(b/=C, runtime_error);
//     // solution
//     C = b;
//     EXPECT_NO_THROW(C /= A);
//     EXPECT_EQ(C*A, b);

//     // division double by matrix
//     EXPECT_NO_THROW(C = 2.5/A);
//     EXPECT_EQ(C*A, 2.5*IdMat(4));
//     EXPECT_NO_THROW(C = 2.5/A(ALL, {1,3}));
//     EXPECT_EQ(C*A(ALL, {1,3}), 2.5*IdMat(3));
//     EXPECT_NO_THROW(C = 2.5/A({1,3}, ALL));
//     EXPECT_EQ(A({1,3}, ALL)*C, 2.5*IdMat(3));
// }

// TEST(Matrix, solve_ls){
//     Matrix A,b,C;
//     A = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,2,1,-1}
//     );
//     b = Matrix(4,1, {2,4,1,3});

//     // shape doesn't match
//     EXPECT_THROW(Matrix<T>::solve_ls(A, b({1,3},0)), invalid_argument);
//     // underdetermined
//     C = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,0,0,0}
//     );
//     EXPECT_THROW(Matrix<T>::solve_ls(C, b), runtime_error);

//     // solution
//     EXPECT_NO_THROW(C = Matrix<T>::solve_ls(A,b));
//     EXPECT_EQ(A * C, b);

//     // overdetermined
//     A = Matrix (6,4,
//         {-1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,2,1,-1,
//          7,-3,0,3,
//          -2,5,1,4}
//     );
//     b = Matrix(6,1, {2,4,1,3,-3,5});
//     EXPECT_NO_THROW(C = Matrix<T>::solve_ls(A,b));
//     Matrix<T>::set_double_precision(8);
//     EXPECT_EQ(C, 
//         Matrix(4,1,{0.34582912, 
//                     1.42291125, 
//                     -0.110443, 
//                     -0.456683})
//     );
//     Matrix<T>::set_double_precision();
// }

// TEST(Matrix, matrix_r_divide){
//     Matrix A,b,C;
//     A = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,2,1,-1}
//     );
//     b = Matrix(1,4, {2,4,1,3});

//     // shape doesn't match
//     EXPECT_THROW(Matrix<T>::matrix_r_divide(b(0,{1,3}), A), invalid_argument);
//     // underdetermined
//     C = Matrix (4,4,
//         {1,3,4,2,
//          0,2,1,-2,
//          2,1,-3,2,
//          0,0,0,0}
//     );
//     EXPECT_THROW(Matrix<T>::matrix_r_divide(b, C), runtime_error);
    
//     // solution
//     EXPECT_NO_THROW(C = Matrix<T>::matrix_r_divide(b, A));
//     EXPECT_EQ(C*A, b);

//     // overdetermined
//     A = Matrix (6,4,
//         {-1,3,4,2,
//         0,2,1,-2,
//         2,1,-3,2,
//         0,2,1,-1,
//         7,-3,0,3,
//         -2,5,1,4}
//     ).t();
//     b = Matrix(1,6, {2,4,1,3,-3,5});
//     EXPECT_NO_THROW(C = Matrix<T>::matrix_r_divide(b,A));
//     Matrix<T>::set_double_precision(8);
//     EXPECT_EQ(C, 
//         Matrix(1,4,{0.34582912, 
//                     1.42291125, 
//                     -0.110443, 
//                     -0.456683})
//     );
//     Matrix<T>::set_double_precision();

// }



// TEST(Matrix, eigenvalues){
//     Matrix m1, D, V;

//     m1 = Matrix (4,4,
//         {7.5231e-01,   8.7419e-01,   3.6122e-01,   4.6593e-01,
//          6.4349e-01,   2.9453e-01,   4.3203e-01,   8.8371e-03,
//          1.4175e-01,   8.3325e-01,   6.4892e-01,   8.3927e-02,
//          8.1433e-01,   9.5796e-01,   9.0255e-01,   1.0307e-01}
//     );

//     // matrix not square
//     EXPECT_THROW(m1({1, m1.r()-1}, ALL).eigen_QR(D, V), invalid_argument);

//     // eigen problem solution
//     EXPECT_NO_THROW(m1.eigen_QR(D, V));

//     Matrix<T>::set_double_precision(6);
//     EXPECT_EQ( D,
//         Matrix(4, 1, {1.90595741, -0.44608142,  0.35662657, -0.01767256}));
//     EXPECT_EQ( V,
//         -Matrix(4, 4, {-0.61559328, -0.4662597 , -0.5442746 , -0.13613666,
//                        -0.33941002,  0.68239286, -0.20427391, -0.46868028,
//                        -0.33624531, -0.48129845,  0.77632961,  0.52722921,
//                        -0.62672549,  0.29205081,  0.24361787,  0.69558246}));
//     EXPECT_EQ( V * diag(D) * V.inv(), m1);

//     /*

//     // testing for complex eigenvalues
//     m1 = RandMat(4,4);
//     cout << m1 << endl;

//     EXPECT_NO_THROW(m1.eigen_QR(D, V, 100000));
//     cout << "QR without shift:" << endl;
//     cout << D << endl;
//     cout << V << endl;

//     EXPECT_NO_THROW(m1.eigen_QR_shift(D, V, 100000));
//     cout << "QR with shift:" << endl;
//     cout << D << endl;
//     cout << V << endl;

//     */

//     Matrix<T>::set_double_precision();

// }

// TEST(Matrix, implicit_double_QR_step){
//     Matrix m1 = RandMat(6,6);
//     Matrix Q,A;
//     m1.hessenberg_dec(Q,A);

//     for(uint i=0; i<5; i++){
//         cout << "step n# " << i << endl;
//         cout << "starting mat: " << A << endl;
//         A = A.implicit_double_QR_step();
//     }
//     cout << "result: " << A << endl;
    
// } 