#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "../matrix_operations.h"

struct Point {
	long double x;
	long double y;
};

std::istream& operator >> (std::istream& is, Point& p) {
	is >> p.x >> p.y;
	return is;
}

struct SplinePolynomial {
	std::vector<long double> coefficients;
	std::pair<long double, long double> boundaries;
};

long double calculateH(const std::vector<Point>& points, size_t i) {
	return points[i].x - points[i - 1].x;
}

std::tuple<Matrix<long double>, Matrix<long double>> generateTridiagonalMatrixEquation(const std::vector<Point>& points) {
	size_t pointsQuantity = points.size();
	size_t equationsNumber = pointsQuantity - 2;
	Matrix<long double> tridiagonalMatrix(equationsNumber, 3);
	Matrix<long double> rightSide(equationsNumber, 1);
	tridiagonalMatrix[0][0] = 0;
	tridiagonalMatrix[0][1] = 2 * (calculateH(points, 1) + calculateH(points, 2));
	tridiagonalMatrix[0][2] = calculateH(points, 2);
	rightSide[0][0] = 3 * ((points[2].y - points[1].y) / calculateH(points, 2) - (points[1].y - points[0].y) / calculateH(points, 2));
	for (size_t row = 1; row + 1 < equationsNumber; ++row) {
		long double hLeft = calculateH(points, row + 1);
		long double hRight = calculateH(points, row + 2);
		tridiagonalMatrix[row][0] = hLeft;
		tridiagonalMatrix[row][1] = 2 * (hLeft + hRight);
		tridiagonalMatrix[row][2] = hRight;
		rightSide[row][0] = 3 * ((points[row + 2].y - points[row + 1].y) / hRight - (points[row + 1].y - points[row].y) / hLeft);
	}
	tridiagonalMatrix[pointsQuantity - 3][0] = calculateH(points, pointsQuantity - 2);
	tridiagonalMatrix[pointsQuantity - 3][1] = 2 * (calculateH(points, pointsQuantity - 2) 
						   + calculateH(points, pointsQuantity - 1));
	rightSide[equationsNumber - 1][0] = 3 * ((points[pointsQuantity - 1].y - points[pointsQuantity - 2].y) / calculateH(points, pointsQuantity - 1)
					    - (points[pointsQuantity - 2].y - points[pointsQuantity - 3].y) / calculateH(points, pointsQuantity - 2));
	return std::make_tuple(tridiagonalMatrix, rightSide);
}

Matrix<long double> findSplinePolynomial(const std::vector<Point>& points) {
	size_t pointsQuantity = points.size();
	Matrix<long double> splineCoefficients(pointsQuantity - 1, 4);
	std::tuple<Matrix<long double>, Matrix<long double>> tridiagonalEquation = generateTridiagonalMatrixEquation(points);
	Matrix<long double> equationSolution = solveTridiagonalMatrixSystemOfLineralEquations(std::get<0>(tridiagonalEquation), 
											      std::get<1>(tridiagonalEquation));
	for (size_t row = 0; row < pointsQuantity - 1; ++row) {
		splineCoefficients[row][0] = points[row].y;
		splineCoefficients[row][2] = (row != 0 ? equationSolution[row - 1][0] : 0);
	}
	for (size_t row = 0; row < pointsQuantity - 1; ++row) {
		if (row + 1 != pointsQuantity - 1) {
			splineCoefficients[row][1] = (points[row + 1].y - points[row].y) / calculateH(points, row + 1) 
						     - (calculateH(points, row + 1) * 
						     (splineCoefficients[row + 1][2] + 2 * splineCoefficients[row][2])) / 3;
			splineCoefficients[row][3] = (splineCoefficients[row + 1][2] - splineCoefficients[row][2]) / (3 * calculateH(points, row + 1));
		} else {
			splineCoefficients[row][1] = (points[row + 1].y - points[row].y) / calculateH(points, row + 1) 
					   	     - 2 * calculateH(points, row + 1) * splineCoefficients[row][2] / 3;
			splineCoefficients[row][3] = -splineCoefficients[row][2] / (3 * calculateH(points, row + 1));
		}
	}
	return splineCoefficients;
}

long double calculateSplinePointValue(const Matrix<long double>& coefficients, const std::vector<Point>& points, long double x) {
	size_t pointsQuantity = points.size();
	std::pair<long double, long double> pointIndexes = {-1, -1};
	for (size_t pointIndex = 1; pointIndex + 1 < pointsQuantity; ++pointIndex) {
		if (x >= points[pointIndex - 1].x && x <= points[pointIndex].x) {
			pointIndexes.first = (int)pointIndex - 1;
			pointIndexes.second = pointIndex;
			break;
		}
	}
	if (pointIndexes.first == -1) {
		throw std::invalid_argument("There are no such segments");
	}
	size_t desiredIndex = pointIndexes.first;
	std::cout << desiredIndex << std::endl;
	return coefficients[desiredIndex][0] 
		+ coefficients[desiredIndex][1] * (x - points[desiredIndex].x) 
		+ coefficients[desiredIndex][2] * std::pow(x - points[desiredIndex].x, 2) 
		+ coefficients[desiredIndex][3] * std::pow(x - points[desiredIndex].x, 3);
	
}

int main() {
	size_t pointsQuantity;
	std::cout << "Enter number of points: " << std::endl;
	std::cin >> pointsQuantity;
	std::vector<Point> points(pointsQuantity);
	std::cout << "Enter points: " << std::endl;
	for (size_t pointIndex = 0; pointIndex < pointsQuantity; ++pointIndex) {
		std::cin >> points[pointIndex];
	}
	std::cout << "Spline coefficients: " << std::endl;
	Matrix<long double> coefficients = findSplinePolynomial(points);
	std::cout << coefficients << std::endl;
	std::cout << "Enter X*" << std::endl;
	long double checkVariable;
	std::cin >> checkVariable;
	std::cout << "Spline value: " << std::endl;
	std::cout << calculateSplinePointValue(coefficients, points, checkVariable) << std::endl;
	return 0;
}
