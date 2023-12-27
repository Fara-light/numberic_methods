#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "../../matrix.h"
#include "../../matrix_operations.h"

std::ostream& operator << (std::ostream& os, const std::complex<long double>& val) {
	long double real = std::real(val);
	long double imag = std::imag(val);
	os << real;
	if (imag) {
		std::cout << (imag < 0 ? " - " : " + ") << std::fabs(imag) << "i";
	}
	return os;
}

long double getVectorEuclideanNorm(const Matrix<long double>& v) {
  return std::sqrt((transpose(v) * v)[0][0]);
}

Matrix<long double> getNormalizedVector(size_t row_number) {
  Matrix<long double> v = Matrix<long double>(row_number, 1);
  for (size_t i = 0; i < row_number; ++i) {
    v[i][0] = 1;
  }
  return v / getVectorEuclideanNorm(v);
}

void fillMatrixColumnWithVector(Matrix<long double>& matrix, Matrix<long double>& vector, size_t column) {
  for (size_t row = 0; row < matrix.getDimention().rows_number; ++row) {
    matrix[row][column] = vector[row][0];
  }
}

std::pair<Matrix<long double>, Matrix<long double>> executeLanczosAlgorithm(const Matrix<long double>& A) {
  Dimention dimention = A.getDimention();
  if (dimention.rows_number != dimention.columns_number) {
    throw std::invalid_argument("Recieved not a square matrix");
  }
  Matrix<long double> T(dimention);
  Matrix<long double> V(dimention);
  Matrix<long double> v = getNormalizedVector(dimention.rows_number);
  fillMatrixColumnWithVector(V, v, 0);
  Matrix<long double> raw_w = A * v;
  long double alpha = (transpose(v) * raw_w)[0][0];
  T[0][0] = alpha;
  Matrix<long double> w = raw_w - alpha * v;
  for (size_t i = 1; i < dimention.rows_number; ++i) {
    Matrix<long double> prev_w = w;
    Matrix<long double> prev_v = v;
    long double beta = getVectorEuclideanNorm(w);
    v = prev_w / beta;
    fillMatrixColumnWithVector(V, v, i);
    raw_w = A * v;
    alpha = (transpose(v) * raw_w)[0][0];
    w = raw_w - alpha * v - beta * prev_v;
    T[i][i] = alpha;
    T[i - 1][i] = beta;
    T[i][i - 1] = beta;
  }
  return std::make_pair(T, V);
}

Matrix<long double> normalizeColumn(const Matrix<long double>& matrix, size_t column) {
  Dimention dimention = matrix.getDimention();
  Matrix<long double> normalizedColumn(dimention.rows_number, 1);
  long double norm = 0;
  for (size_t row = 0; row < dimention.rows_number; ++row) {
    norm += matrix[row][column] * matrix[row][column];
  }
  norm = std::sqrt(norm);
  for (size_t row = 0; row < dimention.rows_number; ++row) {
    normalizedColumn[row][0] = matrix[row][column] / norm;
  }
  return normalizedColumn;
}

Matrix<long double> inverseIteration(const Matrix<long double>& matrix, const std::complex<long double>& targetEigenvalue, long double EPS) {
  Dimention dimention = matrix.getDimention();
  Matrix<long double> shiftedMatrix = matrix - targetEigenvalue.real() * getIdentityMatrix<long double>(dimention);
  while (!qrEndCondition(shiftedMatrix, EPS)) {
    shiftedMatrix = executeIterationQR(shiftedMatrix);
  }
  return normalizeColumn(shiftedMatrix, dimention.columns_number - 1);
}

std::vector<std::pair<std::complex<long double>, Matrix<long double>>> findMatrixEigenPairsQR(const Matrix<long double>& matrix, long double EPS) {
  Dimention dimention = matrix.getDimention();
  Matrix<long double> A = matrix;
  std::vector<std::pair<std::complex<long double>, Matrix<long double>>> eigenPairs;
  while (!qrEndCondition(A, EPS)) {
    A = executeIterationQR(A);
  }
  for (size_t column = 0; column < dimention.columns_number; ++column) {
    if (calculateColumnNormUnderRow(A, column, column) > EPS) {
      auto blockEigenValues = findEigenValuesOf2x2Matrix(copyBlock2x2(A, column, column));
      for (const auto& eigenvalue : blockEigenValues) {
        auto eigenvector = inverseIteration(A, eigenvalue, EPS);
        eigenPairs.push_back({eigenvalue, eigenvector});
      }
      ++column;
    } else {
      std::complex<long double> eigenvalue(A[column][column], 0);
      Matrix<long double> eigenvector = normalizeColumn(A, column);
      eigenPairs.push_back({eigenvalue, eigenvector});
    }
  }
  return eigenPairs;
}

Matrix<long double> getMatrixColumn(const Matrix<long double>& A, size_t column) {
  Dimention dimention = A.getDimention();
  Matrix<long double> v(dimention.rows_number, 1);
  for (size_t row = 0; row < dimention.rows_number; ++row) {
    v[row][0] = A[row][column];
  }
  return v;
}

void orthogonolizeMatrixColumns(Matrix<long double>& A) {
  Dimention dimention = A.getDimention();
  for (size_t column = 1; column < dimention.columns_number; ++column) {
    Matrix<long double> currentColumn = getMatrixColumn(A, column);
    Matrix<long double> deductable(currentColumn.getDimention());
    for (size_t j = 0; j < column; ++j) {
      deductable = deductable + (transpose(currentColumn) * getMatrixColumn(A, j))[0][0] * getMatrixColumn(A, j);
    }
    currentColumn = currentColumn + deductable;
    fillMatrixColumnWithVector(A, currentColumn, column);
  }
}

int main() {
  std::cout << "Enter n from matrix n * n: " << std::endl;
  size_t n;
  std::cin >> n;
  Matrix<long double> A(n, n);
  std::cout << "enter matrix: \n";
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      std::cin >> A[i][j];
    }
  }
  std::cout << std::endl;
  Matrix<long double> T(n, n);
  Matrix<long double> V(n, n);
  std::tie(T, V) = executeLanczosAlgorithm(A);
  std::cout << "Lanczos algorithm result: " << std::endl;
  std::cout << "T: " << std::endl << T << std::endl;
  std::cout << "V: " << std::endl << V << std::endl;
  orthogonolizeMatrixColumns(V);
  orthogonolizeMatrixColumns(V);
  std::cout << "V after orthogonolization: " << std::endl << V << std::endl;
  auto eigenPairs = findMatrixEigenPairsQR(T, 0.0000000001);
  for (const auto& eigenPair: eigenPairs) {
    std::cout << "Eigen value: " << eigenPair.first << std::endl;
    std::cout << "Eigen vector: " << std::endl;
    std::cout << V * eigenPair.second << std::endl;
  }
}
