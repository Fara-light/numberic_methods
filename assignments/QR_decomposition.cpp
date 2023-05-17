#include <iostream>
#include "../matrix.h"
#include "../matrix_operations.h"

template <typename T>
std::ostream& operator << (std::ostream& os, const std::complex<T>& val) {
	T real = std::real(val);
	T imag = std::imag(val);
	os << real;
	if (imag) {
		std::cout << (imag < 0 ? " - " : " + ") << std::fabs(imag) << "i";
	}
	return os;
}

int main() {
	std::cout << "Enter matrix order: ";
	size_t n;
	std::cin >> n;
	long double EPS;
	std::cout << "Enter EPS: ";
	std::cin >> EPS;
	Matrix<long double> matrix(n, n);
	std::cout << "Enter matrix: " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		for (size_t column = 0; column < n; ++column) {
			std::cin >> matrix[row][column];
		}
	}	
	auto ans = findMatrixEigenValuesQR(matrix, EPS);
	std::cout << "Eigen values: " << std::endl;
	for (const auto& eigenValue: ans) {
		std::cout << eigenValue << std::endl;
	}
	std::cout << std::endl;
	return 0;
}
