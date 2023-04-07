#include <iostream>
#include </home/fara/numeric_methods/matrix.h>
#include </home/fara/numeric_methods/matrix_operations.h>

int main() {
	std::cout << "Enter matrix order: ";
	size_t n;
	std::cin >> n;
	Matrix<long double> matrix(n, n);
	std::cout << "Enter matrix: " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		for (size_t column = 0; column < n; ++column) {
			std::cin >> matrix[row][column];
		}
	}
	auto pairLU = decompositionLU(matrix);
	std::cout << "LU-decomposition executed" << std::endl;
	std::cout << "L-matrix: " << std::endl;
	std::cout << pairLU.first << std::endl;
	std::cout << "U-matrix: " << std::endl;
	std::cout << pairLU.second << std::endl;
	Matrix<long double> b(n, 1);
	std::cout << "Enter right side: " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		std::cin >> b[row][0];
	}
	auto y = solveSystemOfLineralEquationsLU(pairLU, b);
	std::cout << "Answer: " << std::endl;
	std::cout << y << std::endl;
	Matrix<long double> inverseMatrix = findInverseMatrix(matrix);
	std::cout << "Inverse matrix: " << std::endl;
	std::cout << inverseMatrix << std::endl;
	std::cout << "Matrix * InverseMatrix: " << std::endl;
	std::cout << matrix * inverseMatrix << std::endl;
	std::cout << "Determinant of matrix: " << std::endl;
	std::cout << "LU: " << calculateDeterminantLU(pairLU) << " Ordinary: " << calculateDeterminant(matrix) << std::endl;
	return 0;
}

