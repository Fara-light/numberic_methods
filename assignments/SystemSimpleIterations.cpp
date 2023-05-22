#include <iostream>
#include <cmath>
#include "../matrix.h"
#include "../matrix_operations.h"

long double functionFirst(const Matrix<long double>& X) {
	return cos(X[1][0]) + 3;
}

long double functionSecond(const Matrix<long double>& X) {
	return sin(X[0][0]) + 3;
}

std::pair<Matrix<long double>, size_t> simpleIterationsMethod(const Matrix<long double (*)(const Matrix<long double>&)>& functionMatrix,
			  	       Matrix<long double> X, long double EPS) {
	size_t iterationsCnt = 0;
	Matrix<long double> prevX = X;
	Matrix<long double> deltaX = X;
	do {
		++iterationsCnt;
		prevX = X;
		X = functionMatrix(prevX);
	} while (calculateVectorNorm(X - prevX) > EPS);
	return std::make_pair(X, iterationsCnt);
}

template<typename T>
Matrix<T> findFirstX(const Matrix<T>& left, const Matrix<T>& right) {
	Matrix<T> result(left.getDimention());
	result[0][0] = left[0][0] + (right[0][0] - left[0][0]) / 2;
	result[1][0] = left[0][0] + (right[1][0] - left[1][0]) / 2;
	return result;
}

int main() {
	Matrix<long double (*)(const Matrix<long double>&)> matrix(2, 1);
	matrix[0][0] = functionFirst;
	matrix[1][0] = functionSecond;
	Matrix<long double (*)(const Matrix<long double>&)> J(2, 2);
	Matrix<long double> left(matrix.getDimention()), right(matrix.getDimention());
	std::cout << "Enter left approximation: " << std::endl;
	for (size_t row = 0; row < matrix.getDimention().rows_number; ++row) {
		std::cin >> left[row][0];
	}
	std::cout << "Enter right approximation: " << std::endl;
	for (size_t row = 0; row < matrix.getDimention().rows_number; ++row) {
		std::cin >> right[row][0];
	}
	std::cout << "Enter EPS: ";
	long double EPS;
	std::cin >> EPS;
	auto X = findFirstX(left, right);
	auto ans = simpleIterationsMethod(matrix, X, EPS);
	std::cout << "Simple iterations method finished in " << ans.second << " iterations" << std::endl;
	std::cout << ans.first;
	return 0;
}
