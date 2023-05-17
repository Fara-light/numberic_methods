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
	Matrix<long double> b(n, 1);
	std::cout << "Enter right side: " << std::endl;
	for (size_t row = 0; row < n; ++row) {
		std::cin >> b[row][0];
	}
	auto pairIterationsCntSolution = solveSystemOfLineralEquationsSimpleIterations(matrix, b, EPS);
	Matrix<long double> x = pairIterationsCntSolution.second;
	size_t iterationsCnt = pairIterationsCntSolution.first;
	std::cout << "Simple iterations method finished in " << iterationsCnt << " iterations" << std::endl;
	std::cout << "Answer: " << std::endl << x << std::endl;
	return 0;
}
