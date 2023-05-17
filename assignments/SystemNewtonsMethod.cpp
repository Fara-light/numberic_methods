#include <iostream>
#include <cmath>
#include "../matrix.h"
#include "../matrix_operations.h"

long double functionFirst(const Matrix<long double>& X) {
	return X[0][0] - cos(X[1][0]) - 3;
}

long double functionSecond(const Matrix<long double>& X) {
	return X[1][0] - sin(X[0][0]) - 3;
}

long double dFunctionFirstdx1(const Matrix<long double>& X) {
	return 1;
}

long double dFunctionFirstdx2(const Matrix<long double>& X) {
	return sin(X[1][0]);
}

long double dFunctionSeconddx1(const Matrix<long double>& X) {
	return -cos(X[0][0]);
}

long double dFunctionSeconddx2(const Matrix<long double>& X) {
	return 1;
}

template <typename T>
Matrix<T> findNewtonsA1(const Matrix<T>& fMatrix, const Matrix<T>& J) {
	Matrix<T> resultingMatrix = J;
	resultingMatrix[0][0] = fMatrix[0][0];
	resultingMatrix[1][0] = fMatrix[1][0];
	return resultingMatrix;
}

template <typename T>
Matrix<T> findNewtonsA2(const Matrix<T>& fMatrix, const Matrix<T>& J) {
	Matrix<T> resultingMatrix = J;
	resultingMatrix[0][1] = fMatrix[0][0];
	resultingMatrix[1][1] = fMatrix[1][0];
	return resultingMatrix;
}

std::pair<Matrix<long double>, size_t> newtonsMethod(const Matrix<long double (*)(const Matrix<long double>&)>& functionMatrix,
	       		  const Matrix<long double (*) (const Matrix<long double>&)>& J, 
			  Matrix<long double> X, 
			  long double EPS) {
	size_t iterationsCnt = 0;
	Matrix<long double> prevX = X;
	Matrix<long double> deltaX = X;
	do {
		++iterationsCnt;
		prevX = X;
		long double detJ = calculateDeterminant(J(X));
		Matrix<long double (*)(const Matrix<long double>&)> A1 = findNewtonsA1(functionMatrix, J);
		Matrix<long double (*)(const Matrix<long double>&)> A2 = findNewtonsA2(functionMatrix, J);
		long double detA1 = calculateDeterminant(A1(X));
		long double detA2 = calculateDeterminant(A2(X));
		deltaX[0][0] = detA1 / detJ;
		deltaX[1][0] = detA2 / detJ;
		X = prevX - deltaX;
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
	J[0][0] = dFunctionFirstdx1;
	J[0][1] = dFunctionFirstdx2;
	J[1][0] = dFunctionSeconddx1;
	J[1][1] = dFunctionSeconddx2;
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
	auto ans = newtonsMethod(matrix, J, X, EPS);
	std::cout << "Newton's method finished in " << ans.second << " iterations" << std::endl;
	std::cout << ans.first;
	return 0;
}
