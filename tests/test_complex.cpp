#include <gtest/gtest.h>
#include "complex/complex.hpp"

using namespace MA;

using std::cout, std::endl, std::string;
using std::invalid_argument, std::runtime_error, std::out_of_range;

TEST(ComplexTest, Construction) {
    // Test default constructor
    MA::Complex default_complex;
    EXPECT_EQ(default_complex.real(), 0.0);
    EXPECT_EQ(default_complex.imag(), 0.0);

    // Test constructor with real part
    MA::Complex real_part_complex(3.14);
    EXPECT_EQ(real_part_complex.real(), 3.14);
    EXPECT_EQ(real_part_complex.imag(), 0.0);

    // Test constructor with pair
    std::pair<double, double> complex_pair = std::make_pair(2.0, -1.0);
    MA::Complex pair_complex(complex_pair);
    EXPECT_EQ(pair_complex.real(), 2.0);
    EXPECT_EQ(pair_complex.imag(), -1.0);
}

TEST(ComplexTest, Assignment) {
    MA::Complex complex1(2.0, 3.0);
    MA::Complex complex2;

    // Test assignment operator with Complex
    complex2 = complex1;
    EXPECT_EQ(complex2, complex1);

    // Test assignment operator with double
    complex2 = 5.0;
    EXPECT_EQ(complex2.real(), 5.0);
    EXPECT_EQ(complex2.imag(), 0.0);

    // Test assignment operator with Complex pointer
    MA::Complex* complex_ptr = new MA::Complex(1.0, -2.0);
    complex2 = complex_ptr;
    EXPECT_EQ(complex2, *complex_ptr);

    delete complex_ptr;
}

TEST(ComplexTest, GetterSetter) {
    MA::Complex complex(2.0, 3.0);

    // Test real getter and setter
    EXPECT_EQ(complex.real(), 2.0);
    complex.real(5.0);
    EXPECT_EQ(complex.real(), 5.0);
    complex.real() = 7.0;
    EXPECT_EQ(complex.real(), 7.0);

    // Test imag getter and setter
    EXPECT_EQ(complex.imag(), 3.0);
    complex.imag(-1.0);
    EXPECT_EQ(complex.imag(), -1.0);
    complex.imag() = 3.0;
    EXPECT_EQ(complex.imag(), 3.0);
}

TEST(ComplexTest, Magnitude) {
    MA::Complex complex(3.0, 4.0);
    EXPECT_EQ(complex.magnitude(), 5.0);

    // Test magnitude of complex with only real part
    MA::Complex real_part_complex(7.0);
    EXPECT_EQ(real_part_complex.magnitude(), 7.0);

    // Test magnitude of complex with only imaginary part
    MA::Complex imag_part_complex(0.0, -2.5);
    EXPECT_EQ(imag_part_complex.magnitude(), 2.5);
}

TEST(ComplexTest, Conjugate) {
    MA::Complex complex(2.0, -3.0);
    EXPECT_EQ(complex.conj(), Complex(2.0, 3.0));
}

TEST(ComplexTest, Equality) {
    MA::Complex complex1(2.0, 3.0);
    MA::Complex complex2(2.0, 3.0);
    MA::Complex complex3(4.0, 5.0);

    EXPECT_TRUE(complex1 == complex2);
    EXPECT_FALSE(complex1 == complex3);
    EXPECT_TRUE(Complex() == 0.0);
}

TEST(ComplexTest, Inequality) {
    MA::Complex complex1(2.0, 3.0);
    MA::Complex complex2(2.0, 3.0);
    MA::Complex complex3(4.0, 5.0);

    EXPECT_FALSE(complex1 != complex2);
    EXPECT_TRUE(complex1 != complex3);
    EXPECT_FALSE(Complex() != 0.0);
}

TEST(ComplexTest, Addition) {
    MA::Complex complex1(2.0, 3.0);
    MA::Complex complex2(4.0, 5.0);

    // Test addition with a Complex
    EXPECT_EQ(complex1 + complex2, Complex(6.0, 8.0));

    // Test addition with a double
    EXPECT_EQ(complex1 + 2.5, Complex(4.5, 3.0));

    // Test self addition with a Complex
    complex2 += complex1;
    EXPECT_EQ(complex2, Complex(6.0, 8.0));

    // Test self addition with a double
    complex1 += 2.5;
    EXPECT_EQ(complex1, Complex(4.5, 3.0));
}

TEST(ComplexTest, Subtraction) {
    MA::Complex complex1(5.0, 8.0);
    MA::Complex complex2(2.0, 3.0);

    // Test subtraction with a Complex
    EXPECT_EQ(complex1 - complex2, Complex(3.0, 5.0));

    // Test unary subtraction
    EXPECT_EQ(-complex1, Complex(-5.0, -8.0));

    // Test subtraction with a double
    EXPECT_EQ(complex1 - 2.5, Complex(2.5, 8.0));

    // Test self subtraction with a Complex 
    complex1 -= complex2;
    EXPECT_EQ(complex1, Complex(3.0, 5.0));

    // Test self subtraction with a double 
    complex2 -= 3;
    EXPECT_EQ(complex2, Complex(-1.0, 3.0));
}

TEST(ComplexTest, Multiplication) {
    MA::Complex complex1(2.0, 3.0);
    MA::Complex complex2(4.0, 5.0);

    // Test multiplication with a Complex 
    EXPECT_EQ(complex1 * complex2, Complex(-7.0, 22.0));

    // Test multiplication with a double
    EXPECT_EQ(complex1 * 2.5, Complex(5.0, 7.5));

    // Test self multiplication with a Complex
    complex2 *= complex1;
    EXPECT_EQ(complex2, Complex(-7.0, 22.0));

    // Test self multiplication with a double
    complex1 *= 2.5;
    EXPECT_EQ(complex1, Complex(5.0, 7.5));
}

TEST(ComplexTest, Division) {
    MA::Complex complex1(10.0, 5.0);
    MA::Complex complex2(2.0, 3.0);

    // Test division by a Complex number
    EXPECT_EQ(complex1 / complex2, Complex(35.0/13.0, -20.0/13.0));

    // Test division with a double
    EXPECT_EQ(complex1 / 2.5, Complex(4.0, 2.0));

    // Test division from a double
    EXPECT_EQ(2.5 / complex1, Complex(0.2, -0.1));

    // Test self division by a Complex number
    complex1 /= complex2;
    EXPECT_EQ(complex1, Complex(35.0/13.0, -20.0/13.0));

    // Test self division with a double
    complex2 /= 2.5;
    EXPECT_EQ(complex2, Complex(4.0/5.0, 6.0/5.0));

    // Test division by zero
    EXPECT_THROW(complex1 / Complex(), runtime_error);
    EXPECT_THROW(complex1 / 0, runtime_error);
    EXPECT_THROW(2.5 / Complex(), runtime_error);
    EXPECT_THROW(complex1 /= Complex(), runtime_error);
    EXPECT_THROW(complex1 /= 0.0, runtime_error);
}

TEST(ComplexTest, OutputOperator) {
    // Test output stream operator for Complex
    MA::Complex complex1(2.0, 3.0);
    ::testing::internal::CaptureStdout();
    std::cout << complex1;
    std::string output1 = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(output1, "2+3j");

    MA::Complex complex2(-4.0, 5.0);
    ::testing::internal::CaptureStdout();
    std::cout << complex2;
    std::string output2 = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(output2, "-4+5j");

    // Test output stream operator for Complex pointer
    MA::Complex* complexPtr1 = new MA::Complex(1.5, -2.5);
    ::testing::internal::CaptureStdout();
    std::cout << complexPtr1;
    std::string outputPtr1 = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(outputPtr1, "1.5-2.5j");

    MA::Complex* complexPtr2 = new MA::Complex(-3.0, 0.0);
    ::testing::internal::CaptureStdout();
    std::cout << complexPtr2;
    std::string outputPtr2 = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ(outputPtr2, "-3");

    delete complexPtr1;
    delete complexPtr2;
}