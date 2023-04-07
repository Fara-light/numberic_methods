#include <iostream>
#include </home/fara/numeric_methods/matrix.h>
#include </home/fara/numeric_methods/matrix_operations.h>

int main() {
	std::cout << "Enter matrix order: ";
	size_t n;
	std::cin >> n;
	Matrix<long double> matrix(n, 3);
	std::cout << "Enter diagonals (a, b, c): " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		std::cin >> matrix[row][0] >> matrix[row][1] >> matrix[row][2];
	}
	Matrix<long double> b(n, 1);
	std::cout << "Enter right side: " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		std::cin >> b[row][0];
	}
	Matrix<long double> x = solveTridiagonalMatrixSystemOfLineralEquations(matrix, b);
	std::cout << "Answer: " << std::endl;
	std::cout << x << std::endl;
	return 0;
}
