#ifndef MATRIX_OPERATIONS
#define MATRIX_OPERATIONS

#include </home/fara/numeric_methods/matrix.h>

template <typename T>
Matrix<T> getAdditionalMinor(const Matrix<T>& matrix, size_t rowToRemove, size_t columnToRemove) {
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
	Matrix<T> additionalMinor(dimention.rows_number - 1, dimention.columns_number - 1);
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
			additionalMinor[row - rowDiff][column - columnDiff] = matrix[row][column];
		}
		columnDiff = 0;
	}
	return additionalMinor;
}

#endif
