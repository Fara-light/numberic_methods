#ifndef H_MATRIX_OPERATIONS
#define H_MATRIX_OPERATIONS

#include "/home/fara/numeric_methods/matrix.h"
#include <stdexcept>

template <typename T>
size_t findNonZeroRow(const Matrix<T>& matrix, size_t iteration) {
	// finds first row from from above which has non-zero element on iteration column
	size_t swap_row = iteration + 1;
	size_t dimention = matrix.getDimention().rows_number;
	while (matrix[swap_row][iteration] == 0 && swap_row < dimention) {
		++swap_row;
	}
	return swap_row;
}

template <typename T>
Matrix<T> getLowerDecompositionMatrix(const Matrix<T>& matrix) {
	// returns L-matrix for LU-decomposition
	Dimention matrix_dimention = matrix.getDimention();
	size_t rows_number = matrix_dimention.rows_number;
	size_t columns_number = matrix_dimention.columns_number;
	if (rows_number != columns_number) {
		throw std::exception("invalid argument: rows_number is not equal to columns_number");
	}
	size_t dimention = rows_number;
	Matrix<T> resulting_matrix(rows_number, columns_number);
	for (size_t iteration = 1; iteration < dimention; ++iteration) {
		if (matrix[iteration][iteration] == 0) {
			size_t nonZeroRow = findNonZeroRow(matrix, iteration);
			if (nonZeroRow == dimention) {
				// maybe should throw exception? 
				continue;
			} else {
				swap(matrix[iteration], matrix[nonZeroRow]);
			}
		}
		for (size_t column_number = iteration + 1; column_number < dimention; ++column_number) {
			
		}
	}
}

template <typename T>
Matrix<T> getUpperDecompositionMatrix(const Matrix<T>& matrix) {
	// returns U-matrix for LU-decomposition
}

#endif
