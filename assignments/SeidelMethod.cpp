#include <iostream>
#include </home/fara/numeric_methods/matrix.h>
#include </home/fara/numeric_methods/matrix_operations.h>

int main() {
	size_t n;
	std::cout << "Enter matrix order: ";
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
	Matrix<long double> b(n, 1);
	std::cout << "Enter right side of equation: " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		std::cin >> b[row][0];
	}
	std::pair<size_t, Matrix<long double>> solution = solveSystemOfLineralEquationsSeidelMethod(matrix, b, EPS);
	std::cout << "Seidel method finished in " << solution.first << " iterations" << std::endl;
	std::cout << "Answer: " << std::endl;
	std::cout << solution.second;
	return 0;
}
