#ifndef MATRIX_OPERATIONS
#define MATRIX_OPERATIONS

#include <cmath>
#include <tuple>
#include <iostream>
#include </home/fara/numeric_methods/matrix.h>

template <typename T>
Matrix<T> getMatrixWithExcludedRowColumn(const Matrix<T>& matrix, size_t rowToRemove, size_t columnToRemove) {
	// calculates additional minor of matrix, removing row and column with numbers rowToRemove and columnToRemove
	Dimention dimention = matrix.getDimention();
	if (rowToRemove >= dimention.rows_number) {
		throw std::invalid_argument("row is out of matrix bounds");
	}
	if (columnToRemove >= dimention.columns_number) {
		throw std::invalid_argument("column is out of matrix bounds");
	}
	if (dimention.rows_number <= 1 || dimention.columns_number <= 1) {
		throw std::logic_error("cant create matrix with no elements");
	}
	Matrix<T> matrixWithExcludedRowColumn(dimention.rows_number - 1, dimention.columns_number - 1);
	size_t rowDiff = 0;
	size_t columnDiff = 0;
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		if (row == rowToRemove) {
			rowDiff = 1;
			continue;
		}
		for (size_t column = 0; column < dimention.columns_number; ++column) {
			if (column == columnToRemove) {
				columnDiff = 1;
				continue;
			}
			matrixWithExcludedRowColumn[row - rowDiff][column - columnDiff] = matrix[row][column];
		}
		columnDiff = 0;
	}
	return matrixWithExcludedRowColumn;
}

template <typename T>
T _calculateDeterminant(const Matrix<T>& matrix) {
	// auxiliary function for calculateDeterminant(const Matrix<T>&) function
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number == 1) {
		return matrix[0][0];
	}
	T determinant = 0;
	int multiplier = 1;
	for (size_t column = 0; column < dimention.columns_number; ++column) {
		determinant += multiplier * matrix[0][column] * _calculateDeterminant(getMatrixWithExcludedRowColumn(matrix, 0, column));
		multiplier *= -1;
	}
	return determinant;
}

template <typename T>
T calculateDeterminant(const Matrix<T>& matrix) {
	// calculates determinant of matrix
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't calculate determinant of non-square matrix");
	}
	return _calculateDeterminant(matrix);
}

template <typename T>
std::pair<Matrix<T>, Matrix<T>> decompositionLU(const Matrix<T>& matrix) {
	// executes LU-decomposition of Matrix<T>
	// returns pair of matrices L and U
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't execute decomposition of non-square matrix");
	}
	Matrix<T> U = matrix;
	Matrix<T> L(dimention);
	for (size_t index = 1; index < dimention.rows_number; ++index) {
		for (size_t column = index - 1; column < dimention.columns_number; ++column) {
			for (size_t row = column; row < dimention.rows_number; ++row) {
				if (U[column][column] == 0) {
					throw std::invalid_argument("can't divide by zero");
				}
				L[row][column] = U[row][column] / U[column][column];
			}
		}
		for (size_t row = index; row < dimention.rows_number; ++row) {
			for (size_t column = index - 1; column < dimention.columns_number; ++column) {
				U[row][column] = U[row][column] - L[row][index - 1]*U[index - 1][column];
			}
		}
	}
	return std::make_pair(L, U);
}

template <typename T>
Matrix<T> executeForwardSubstitution(const Matrix<T>& L, const Matrix<T>& b) {
	// auxiliary function, executes forward substitution
	// left side: L, right side: b
	Dimention dimention = L.getDimention();
	size_t n = dimention.rows_number;
	Matrix<T> y(n, 1);
	for (size_t row = 0; row < n; ++row) {
		T prevSum = 0;
		for (size_t column = 0; column < row; ++column) {
			prevSum += L[row][column] * y[column][0];
		}
		y[row][0] = (b[row][0] - prevSum) / L[row][row];
	}
	return y;
}

template <typename T>
Matrix<T> executeBackSubstitution(const Matrix<T>& U, const Matrix<T>& y) {
	// auxiliary function, executes back substitution
	// left side: U, right side: y
	Dimention dimention = U.getDimention();
	size_t n = dimention.rows_number;
	Matrix<T> x(n, 1);
	// overflow is my best friend
	for (size_t row = n - 1; row < n; --row) {
		T prevSum = 0;
		for (size_t column = n - 1; column > row; --column) {
			prevSum += U[row][column] * x[column][0];
		}
		x[row][0] = (y[row][0] - prevSum) / U[row][row];
	}
	return x;
}

template <typename T>
Matrix<T> solveSystemOfLineralEquationsLU(const std::pair<Matrix<T>, Matrix<T>>& pairLU, const Matrix<T>& b) {
	// solves system of lineral equations Ax = b, there A = LU
	// returns vector x
	Dimention bDimention = b.getDimention();
	Dimention aDimention = pairLU.first.getDimention();
	if (bDimention.rows_number != aDimention.rows_number && bDimention.columns_number != 1) {
		throw std::invalid_argument("wrong right side of equation");
	}
	Matrix<T> y = executeForwardSubstitution(pairLU.first, b);
	Matrix<T> x = executeBackSubstitution(pairLU.second, y);
	return x;
}

template <typename T>
T calculateDeterminantLU(const std::pair<Matrix<T>, Matrix<T>>& pairLU) {
	// calculates determinants of A = L * U
	Matrix<T> U = pairLU.second;
	Dimention dimention = U.getDimention();
	T determinant = 1;
	for (size_t index = 0; index < dimention.rows_number; ++index) {
		determinant *= U[index][index];
	}
	return determinant;
}

template <typename T>
void setMatrixColumn(Matrix<T>& matrix, const Matrix<T>& column, size_t column_index) {
	// auxiliary function, replaces matrix column with column_index with input column 
	Dimention dimention = matrix.getDimention();
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		matrix[row][column_index] = column[row][0];
	}
}

template <typename T>
Matrix<T> findInverseMatrix(const Matrix<T>& matrix) {
	// calculates inverse matrix using LU decomposition
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't calculate inverse matrix of non-square matrix");
	}
	if (calculateDeterminant(matrix) == 0) {
		throw std::invalid_argument("can't calculate inverse matrix of degenerate matrix");
	}
	Matrix<T> inverseMatrix(dimention);
	std::pair<Matrix<T>, Matrix<T>> pairLU = decompositionLU(matrix);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		Matrix<T> singleOneVector(dimention.rows_number, 1);
		singleOneVector[row][0] = 1;
		Matrix<T> anotherColumn = solveSystemOfLineralEquationsLU(pairLU, singleOneVector);
		setMatrixColumn(inverseMatrix, anotherColumn, row);
	}
	return inverseMatrix;
}

template <typename T>
std::vector<std::pair<T, T>> calculateRunThroughCoefficients(const Matrix<T>& matrix, const Matrix<T>& B) {
	Dimention dimention = matrix.getDimention();
	std::vector<std::pair<T, T>> runThroughCoefficients;
	runThroughCoefficients.emplace_back(-matrix[0][2]/matrix[0][1], B[0][0]/matrix[0][1]);
	for (size_t row = 1; row < dimention.rows_number; ++row) {
		T prevP = runThroughCoefficients.back().first;
		T prevQ = runThroughCoefficients.back().second;
		T a = matrix[row][0], b = matrix[row][1], c = matrix[row][2];
		T P = -c/(b + a * prevP);
		T Q = (B[row][0] - a * prevQ)/(b + a * prevP);
		runThroughCoefficients.emplace_back(P, Q);
	}
	return runThroughCoefficients;
}

template <typename T>
Matrix<T> solveTridiagonalMatrixSystemOfLineralEquations(const Matrix<T>& matrix, const Matrix<T>& b) {
	// solves system of lineral equations defined by tridiagonal matrix
	Dimention m_dimention = matrix.getDimention();
	Dimention b_dimention = b.getDimention();
	if (m_dimention.rows_number != b_dimention.rows_number && b_dimention.columns_number != 1) {
		throw std::invalid_argument("ivalid b-matrix argument");
	}
	size_t n = m_dimention.rows_number;
	std::vector<std::pair<T, T>> runThroughCoefficients = calculateRunThroughCoefficients(matrix, b);
	Matrix<T> x(n, 1);
	x[n - 1][0] = runThroughCoefficients[n - 1].second;
	for (size_t row = n - 2; row < n; --row) {
		x[row][0] = runThroughCoefficients[row].first * x[row + 1][0] + runThroughCoefficients[row].second;
	}
	return x;
}

template <typename T>
std::pair<Matrix<T>, Matrix<T>> calculateSimpleIterationsCoefficients(const Matrix<T>& matrix, const Matrix<T>& b) {
	// auxiliary function, calculates matrices alpha and beta for simple iterations method
	Dimention dimention = matrix.getDimention();
	Matrix<T> alpha(dimention);
	for(size_t row = 0; row < dimention.rows_number; ++row) {
		for (size_t column = 0; column < dimention.columns_number; ++column) {
			if (row == column) {
				alpha[row][column] = 0;
				continue;
			}
			alpha[row][column] = -matrix[row][column] / matrix[row][row];
		}
	}
	Matrix<T> beta(dimention.rows_number, 1);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		beta[row][0] = b[row][0] / matrix[row][row];
	}
	return std::make_pair(alpha, beta);
}

template <typename T>
T calculateVectorNorm(const Matrix<T>& v) {
	// calculates norm of vector
	Dimention dimention = v.getDimention();
	if (dimention.columns_number != 1) {
		throw std::invalid_argument("can't calculate norm of matrix");
	}
	T norm = 0;
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		norm += v[row][0] * v[row][0];
	}
	return norm;
}

template <typename T>
T calculateMatrixNorm(const Matrix<T>& m) {
	// calculates norm of matrix
	Dimention dimention = m.getDimention();
	T resultingNorm = -1;
	for (size_t column = 0; column <= dimention.columns_number; ++column) {
		T columnSum = 0;
		for (size_t row = 0; row < dimention.rows_number; ++row) {
			columnSum += m[row][column];
		}
		resultingNorm = std::max(resultingNorm, columnSum);
	}
	return resultingNorm;
}

template <typename T>
std::pair<size_t, Matrix<T>> solveSystemOfLineralEquationsSimpleIterations(const Matrix<T>& matrix, const Matrix<T>& b, T EPS) {
	// solves system of lineral equations numerically, using method of simple iterations
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't solve system of lineral equations with non-square matrix");	
	}
	std::pair<Matrix<T>, Matrix<T>> alphaBetaPair = calculateSimpleIterationsCoefficients(matrix, b);
	Matrix<T> alpha = alphaBetaPair.first;
	Matrix<T> beta = alphaBetaPair.second;
	Matrix<T> previousX(dimention.rows_number, 1);
	Matrix<T> currentX = beta;
	size_t iterationsCnt = 0;
	while (fabs(calculateVectorNorm(previousX - currentX)) > EPS) {
		++iterationsCnt;
		std::swap(previousX, currentX);
		currentX = beta + alpha * previousX;
	}
	return std::make_pair(iterationsCnt, currentX);
}

template <typename T>
struct SeidelCoefficients {
	Matrix<T> matrix;
	Matrix<T> vector;
};

template <typename T>
Matrix<T> findInverseOfLowerTriangularMatrix(const Matrix<T>& LT) {
	// finds inverse matrix for lower triangular matrix
	Dimention dimention = LT.getDimention();
	Matrix<T> inverseMatrix(dimention);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		Matrix<T> singleOneVector(dimention.rows_number, 1);
		singleOneVector[row][0] = 1;
		Matrix<T> anotherColumn = executeForwardSubstitution(LT, singleOneVector);
		setMatrixColumn(inverseMatrix, anotherColumn, row);
	}
	return inverseMatrix;
}

template <typename T>
SeidelCoefficients<T> calculateSeidelMethodCoefficients(const std::pair<Matrix<T>, Matrix<T>> alphaBetaPair) {
	// auxiliary function, calculates coefficients for Seidel method
	Matrix<T> alpha = alphaBetaPair.first;
	Matrix<T> beta = alphaBetaPair.second;
	Dimention dimention = alpha.getDimention();
	Matrix<T> C(dimention);
	for (size_t column = 0; column < dimention.columns_number; ++column) {
		for (size_t row = 0; row <= column; ++row) {
			C[row][column] = alpha[row][column];
		}
	}
	Matrix<T> B(dimention);
	for (size_t column = 0; column < dimention.columns_number; ++column) {
		for (size_t row = column + 1; row < dimention.rows_number; ++row) {
			B[row][column] = alpha[row][column];
		}
	}
	Matrix<T> commonPart = findInverseOfLowerTriangularMatrix(getIdentityMatrix<T>(dimention.rows_number) - B);
	return {commonPart * C, commonPart * beta};
}

template <typename T>
std::pair<size_t, Matrix<T>> solveSystemOfLineralEquationsSeidelMethod(const Matrix<T>& matrix, const Matrix<T>& b, T EPS) {
	// solves system of lineral equations with Seidel Method
	// returns std::pair<size_t, Matrix<T>> - number of iterations, answer
	Dimention dimention = matrix.getDimention();
	std::pair<Matrix<T>, Matrix<T>> alphaBetaPair = calculateSimpleIterationsCoefficients(matrix, b);
	SeidelCoefficients<T> seidelCoefficients = calculateSeidelMethodCoefficients(alphaBetaPair);
	Matrix<T> previousX(dimention.rows_number, 1);
	Matrix<T> currentX = alphaBetaPair.second;
	size_t iterationsCnt = 0;
	while (fabs(calculateVectorNorm(previousX - currentX)) >= EPS) {
		std::swap(currentX, previousX);
		currentX = seidelCoefficients.matrix * previousX + seidelCoefficients.vector;
		++iterationsCnt;
	}
	return std::make_pair(iterationsCnt, previousX);
}

template <typename T>
std::pair<size_t, size_t> findLargestAbsOffDiagonalValuePos(const Matrix<T>& matrix) {
	// returns pair of row and column of the largest off diagonal element of matrix by absolute value
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't handle non-square matrix");
	}
	if (dimention.rows_number == 1) {
		throw std::logic_error("can't find max off-diagonal value when no off-diagonal values");
	}
	T maxValue = fabs(matrix[0][1]);
	size_t maxRow = 0;
	size_t maxColumn = 1;
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		for (size_t column = 0; column < dimention.columns_number; ++column) {
			if (row == column) {
				continue;
			}
			if (std::fabs(matrix[row][column]) > maxValue) {
				maxValue = std::fabs(matrix[row][column]);
				maxRow = row;
				maxColumn = column;
			}
		}
	}
	return std::make_pair(maxRow, maxColumn);
}

template <typename T>
T sumSquareOffDiagonalElements(const Matrix<T>& matrix) {
	// calculates sum of off-diagonal elements of matrix
	Dimention dimention = matrix.getDimention();
	T sum = 0;
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		for (size_t column = 0; column < dimention.columns_number; ++column) {
			if (row == column) {
				continue;
			}
			sum += matrix[row][column] * matrix[row][column];
		}
	}
	return sum;
}

template <typename T>
bool rotationMethodEndCondition(const Matrix<T>& matrix, T EPS) {
	// returns true when end-condition is satisfied
	T sum = sumSquareOffDiagonalElements(matrix);
	return std::sqrt(sum) < EPS;
}



template <typename T>
Matrix<T> findRotationMethodMatrixU(const Matrix<T>& matrix, std::pair<size_t, size_t> maxElementRowColumnPair) {
	// finds matrix U for rotation method
	Dimention dimention = matrix.getDimention();
	Matrix<T> U = getIdentityMatrix<T>(dimention);
	size_t maxElementRow = maxElementRowColumnPair.first;
	size_t maxElementColumn = maxElementRowColumnPair.second;
	T phi = std::atan2(2 * matrix[maxElementRow][maxElementColumn], 
			matrix[maxElementRow][maxElementRow] - matrix[maxElementColumn][maxElementColumn]) / 2;
	U[maxElementRow][maxElementRow] = cos(phi);
	U[maxElementRow][maxElementColumn] = -sin(phi);
	U[maxElementColumn][maxElementRow] = sin(phi);
	U[maxElementColumn][maxElementColumn] = cos(phi);
	return U;
}

template <typename T>
Matrix<T> getMatrixColumn(const Matrix<T>& matrix, size_t columnNumber) {
	// returns matrices column with number columnNumber
	Dimention dimention = matrix.getDimention();
	Matrix<T> resultingColumn(dimention.rows_number, 1);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		resultingColumn[row][0] = matrix[row][columnNumber];
	}
	return resultingColumn;
}

template <typename T>
std::tuple<size_t, Matrix<T>, std::vector<Matrix<T>>> findEigenValuesVectorsRotationMethod(const Matrix<T>& matrix, T EPS) {
	// returns number of iterations, vector of eigen values, std::vector of eigen vectors
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't find eigen values of non-square matrix");
	}
	Matrix<T> ans = matrix;
	Matrix<T> eigenVectorMatrixU = getIdentityMatrix<T>(dimention);
	size_t iterationsCount = 0;
	while (not rotationMethodEndCondition(ans, EPS)) {
		std::pair<size_t, size_t> maxElementRowColumnPair = findLargestAbsOffDiagonalValuePos(ans);
		Matrix<T> U = findRotationMethodMatrixU(ans, maxElementRowColumnPair);
		eigenVectorMatrixU = eigenVectorMatrixU * U;
		ans = transpose(U) * ans * U;
		++iterationsCount;
	}
	Matrix<T> eigenValuesVector(dimention.rows_number, 1);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		eigenValuesVector[row][0] = ans[row][row];
	}
	std::vector<Matrix<T>> eigenVectors;
	for (size_t column = 0; column < dimention.columns_number; ++column) {
		eigenVectors.push_back(getMatrixColumn(eigenVectorMatrixU, column));
	}
	return std::make_tuple(iterationsCount, eigenValuesVector, eigenVectors);
}

#endif
