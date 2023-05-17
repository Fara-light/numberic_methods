#include <iostream>
#include "../matrix.h"
#include "../matrix_operations.h"

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
	std::tuple<size_t, Matrix<long double>, std::vector<Matrix<long double>>> solution = findEigenValuesVectorsRotationMethod(matrix, EPS);
	std::cout << "Rotation method finished in " << std::get<0>(solution) << " iterations" << std::endl;
	std::cout << "Ans: " << std::endl;
	std::cout << std::get<1>(solution);
	std::cout << "Eigen vectors: " << std::endl;
	for (const auto& eigenVector: std::get<2>(solution)) {
		std::cout << eigenVector << std::endl;
	}
	return 0;
}
