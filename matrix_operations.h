#ifndef MATRIX_OPERATIONS
#define MATRIX_OPERATIONS

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
	Dimention dimention = matrix.getDimention();
	if (dimention.rows_number != dimention.columns_number) {
		throw std::invalid_argument("can't execute decomposition of non-square matrix");
	}
	Matrix<T> U = matrix;
	Matrix<T> L(dimention);
	for (size_t column = 0; column < dimention.columns_number; ++column) {
		for (size_t row = 0; row < dimention.rows_number; ++row) {
			L[row][column] = U[row][column] / U[column][column];
		}
	}
	for (size_t index = 1; index < dimention.rows_number; ++index) {
		for (size_t column = index - 1; column < dimention.columns_number; ++column) {
			for (size_t row = column; row < dimention.rows_number; ++row) {
				L[row][column] = U[row][column] / U[column][column];
			}
		}
		for (size_t row = index; row < dimention.rows_number; ++row) {
			for (size_t column = index - 1; column < dimention.columns_number; ++column) {
				U[row][column] = U[row][column] - L[row][index - 1]*U[index - 1][column];
			}
		}
	}
	for (size_t column = 1; column < dimention.columns_number; ++column) {
		for (size_t row = 0; row < column; ++row) {
			L[row][column] = 0;
		}
	}
	return std::make_pair(L, U);
}

#endif
