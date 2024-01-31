#include "linear_algebra_ma/matrices.hpp"

int main(){
    using namespace MA;

    std::cout << std::is_arithmetic_v<std::complex<double>> << std::endl;
    std::cout << std::is_arithmetic_v<double> << std::endl;
    Matrix m(2,2);

    return 0;
}