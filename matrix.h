#ifndef H_MATRIX
#define H_MATRIX

#include <iomanip>
#include <algorithm>
#include <utility>
#include <vector>

struct Dimention {
	Dimention() {}
	Dimention(size_t rows_number, size_t columns_number)
		: rows_number(rows_number), columns_number(columns_number) {}
	size_t rows_number = 0;
	size_t columns_number = 0;
};

bool operator == (const Dimention& left, const Dimention& right) {
	return (left.rows_number == right.rows_number) && (left.columns_number == right.columns_number);
}

bool operator != (const Dimention& left, const Dimention& right) {
	return !(left == right);
}

template <typename T>
class Row {
public:
	Row() {}
	Row (size_t n) {
		row.assign(n, 0);
	}
	Row(size_t n, T value) {
		row.assign(n, value);
	}
	size_t getSize() const {
		return row.size();
	}
	T& operator [] (size_t column_number) {
		return row[column_number];
	}
	T operator [] (size_t column_number) const {
		return row[column_number];
	}
private:
	std::vector<T> row;
};

template <typename T>
Row<T> operator * (const Row<T>& row, T multiplier) {
	int columns_number = row.getSize();
	Row<T> resultingRow(row.getSize());
	for (size_t column_number = 0; column_number < columns_number; ++column_number) {
		resultingRow = row[column_number] * multiplier;
	}
	return row;
}

template <typename T>
Row<T> operator * (T multiplier, const Row<T>& row) {
	return row * multiplier;
}

template <typename T>
Row<T> operator + (const Row<T>& left, const Row<T>& right) {
	size_t left_size = left.getSize();
	size_t right_size = right.getSize();
	if (left_size != right_size) {
		throw std::invalid_argument("recieved rows with different size");
	}
	Row<T> resulting_row(left_size);
	for (size_t column_number = 0; column_number < left_size; ++column_number) {
		resulting_row[column_number] = left[column_number] + right[column_number];
	}
	return resulting_row;
}

template <typename T>
class Column {
public:
	Column() {
	}
	Column(size_t n) {
		column.assign(n, 0);
	}
	Column(size_t n, T value) {
		column.assign(n, value);
	}
	size_t getSize() const {
		return column.size();
	}
	T& operator [] (size_t row_number) {
		return column[row_number];
	}
	T operator [] (size_t row_number) const {
		return column[row_number];
	}
private:
	std::vector<T> column;
};

template <typename T>
class Matrix {	
public:
	Matrix(size_t n, size_t m) {
		matrix = Column<Row<T>>(n, Row<T>(m, 0));
	}
	Matrix(const Dimention& dimention) 
		: Matrix(dimention.rows_number, dimention.columns_number) {}
	Dimention getDimention() const {
		return Dimention(matrix.getSize(), matrix[0].getSize());
	}
	void swapRows(size_t left, size_t right) {
		checkRowBounds(left);
		checkRowBounds(right);
		std::swap(matrix[left], matrix[right]);
	}
	void swapColumns(size_t left, size_t right) {
		checkColumnBounds(left);
		checkColumnBounds(right);
		Dimention dimention = getDimention();
		for (size_t row = 0; row < dimention.rows_number; ++row) {
			std::swap(matrix[row][left], matrix[row][right]);
		}
	}
	Row<T>& operator [] (size_t row_number) {
		checkRowBounds(row_number);
		return matrix[row_number];
	}
	Row<T> operator [] (size_t row_number) const {
		checkRowBounds(row_number);
		return matrix[row_number];
	}
private:
	void checkRowBounds(size_t row) const {
		Dimention dimention = getDimention();
		if (dimention.rows_number <= row) {
			throw std::invalid_argument("Row number of first element is out of bounds");
		}
	}
	void checkColumnBounds(size_t column) const {
		Dimention dimention = getDimention();
		if (dimention.columns_number <= column) {
			throw std::invalid_argument("Column number of first element is out of bounds");
		}
	}
	Column<Row<T>> matrix;
};

template <typename T>
Matrix<T> getIdentityMatrix(size_t n) {
	Matrix<T> identityMatrix(n, n);
	for (size_t index = 0; index < n; ++index) {
		identityMatrix[index][index] = 1;
	}
	return identityMatrix;
}

template <typename T>
std::ostream& operator << (std::ostream& os, Matrix<T> matrix) {
	Dimention matrix_dimention = matrix.getDimention();
	for (size_t row_number = 0; row_number < matrix_dimention.rows_number; ++row_number) {
		for (size_t column_number = 0; column_number < matrix_dimention.columns_number; ++column_number) {
			os << std::setw(7) << matrix[row_number][column_number] << " ";
		}
		os << std::endl;
	}
	return os;
}

template <typename T>
Matrix<T> operator * (const Matrix<T>& matrix, T multiplier) {
	Dimention dimention = matrix.getDimention();
	Matrix<T> resultingMatrix(dimention);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		for (size_t column = 0; column < dimention.columns_number; ++column) {
			resultingMatrix[row][column] = matrix[row][column] * multiplier;
		}
	}
	return resultingMatrix;
}

template <typename T>
Matrix<T> operator * (T multiplier, const Matrix<T>& matrix) {
	return matrix * multiplier;
}

template <typename T>
Matrix<T> operator + (const Matrix<T>& left, const Matrix<T>& right) {
	if (left.getDimention() != right.getDimention()) {
		throw std::invalid_argument("sum of matrixes with different dimentions");
	}
	Dimention dimention = left.getDimention();
	Matrix<T> resulting_matrix(dimention);
	for (size_t row = 0; row < dimention.rows_number; ++row) {
		for (size_t column = 0; column < dimention.columns_number; ++column) {
			resulting_matrix[row][column] = left[row][column] + right[row][column];
		}
	}
	return resulting_matrix;
}

template <typename T>
Matrix<T> operator - (const Matrix<T>& left, const Matrix<T>& right) {
	return left + (right * static_cast<T>(-1));
}

template <typename T>
T calculateElementOnPosition(const Matrix<T>& left, const Matrix<T>& right, size_t left_row, size_t right_column) {
	// calculates sum of corresponding product of elements from left matrix row and right matrix column
	Dimention left_dimention = left.getDimention();
	Dimention right_dimention = right.getDimention();
	T resulting_sum = 0;
	for (size_t iterating_pos = 0; iterating_pos < left_dimention.columns_number; ++iterating_pos) {
		resulting_sum = resulting_sum + left[left_row][iterating_pos] * right[iterating_pos][right_column];
	}
	return resulting_sum;
}

template <typename T>
Matrix<T> operator * (const Matrix<T>& left, const Matrix<T>& right) {
	Dimention left_dimention = left.getDimention();
	Dimention right_dimention = right.getDimention();
	if (left_dimention.columns_number != right_dimention.rows_number) {
		throw std::invalid_argument("multiplication of matrixes with wrong dimentions");
	}
	Matrix<T> resulting_matrix(left_dimention.rows_number, right_dimention.columns_number);
	for (size_t left_row = 0; left_row < left_dimention.rows_number; ++left_row) {
		for (size_t right_column = 0; right_column < right_dimention.columns_number; ++right_column) {
			T matrix_element_value = calculateElementOnPosition(left, right, left_row, right_column);
			resulting_matrix[left_row][right_column] = matrix_element_value;
		}
	}
	return resulting_matrix;
}

#endif
