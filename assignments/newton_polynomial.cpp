#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
using namespace std;

struct Point {
	long double x;
	long double y;
};

std::vector<Point> generatePointsByFunctionAndAbscissaCoordinates(long double (*function)(long double), const std::vector<long double>& abscissaCoordinates) {
	std::vector<Point> points;
	for (size_t pointIndex = 0; pointIndex < abscissaCoordinates.size(); ++pointIndex) {
		long double currentCoordinate = abscissaCoordinates[pointIndex];
		points.push_back({currentCoordinate, function(currentCoordinate)});
	}
	return points;
}

std::vector<long double> findNewtonPolynomialCoefficients(const std::vector<Point>& points) {
	size_t coordinatesQuantity = points.size();
	std::vector<std::vector<long double>> dividedDifferenceMatrix(coordinatesQuantity + 1, std::vector<long double>(coordinatesQuantity, 0));
	for (size_t coordinateIndex = 0; coordinateIndex < coordinatesQuantity; ++coordinateIndex) {
		dividedDifferenceMatrix[0][coordinateIndex] = points[coordinateIndex].x;
		dividedDifferenceMatrix[1][coordinateIndex] = points[coordinateIndex].y;
	}
	std::vector<long double> ans(coordinatesQuantity);
	ans[0] = points[0].y;
	for (size_t dividedDifferenceOrder = 1; dividedDifferenceOrder < coordinatesQuantity; ++dividedDifferenceOrder) {
		size_t row = dividedDifferenceOrder + 1;
		for (size_t coordinateIndex = 0; coordinateIndex < coordinatesQuantity - dividedDifferenceOrder; ++coordinateIndex) {
			long double dividedDifference = (dividedDifferenceMatrix[row - 1][coordinateIndex] - dividedDifferenceMatrix[row - 1][coordinateIndex + 1]) / (points[coordinateIndex].x - points[coordinateIndex + dividedDifferenceOrder].x);
			dividedDifferenceMatrix[row][coordinateIndex] = dividedDifference;
		}
		ans[dividedDifferenceOrder] = dividedDifferenceMatrix[row][0];
	}
	for (const auto& line: dividedDifferenceMatrix) {
		for (const auto& element: line) {
			std::cout << element << " ";
		}
		std::cout << std::endl;
	}
	return ans;
}

long double calculatePolynomialValue(const std::vector<long double>& coefficients, const std::vector<long double>& abscissaCoordinates, long double x) {
	if (coefficients.size() != abscissaCoordinates.size()) {
		throw std::invalid_argument("printPolynom different sizes");
	}
	long double polynomialValue = 0;
	for (size_t coordinateIndex = 0; coordinateIndex < abscissaCoordinates.size(); ++coordinateIndex) {
		long double currentMonoidValue = coefficients[coordinateIndex];
		for (size_t pairedCoordinateIndex = 0; pairedCoordinateIndex < coordinateIndex; ++pairedCoordinateIndex) {
			currentMonoidValue *= (x - abscissaCoordinates[pairedCoordinateIndex]);
		}
		polynomialValue += currentMonoidValue;
	}
	return polynomialValue;
}

void printNumberSign(long double number) {
	if (number >= 0) {
		std::cout << " + ";
	} else if (number < 0) {
		std::cout << " - ";
	}
}

void printPolynomial(const std::vector<long double>& coefficients, const std::vector<long double>& abscissaCoordinates) {
	if (coefficients.size() != abscissaCoordinates.size()) {
		throw std::invalid_argument("printPolynom different sizes");
	}
	for (size_t coordinateIndex = 0; coordinateIndex < abscissaCoordinates.size(); ++coordinateIndex) {
		if (coordinateIndex != 0 || coefficients[coordinateIndex] < 0) {
			printNumberSign(coefficients[coordinateIndex]);
		}
		std::cout << std::fabs(coefficients[coordinateIndex]);
		for (size_t pairedCoordinateIndex = 0; pairedCoordinateIndex < coordinateIndex; ++pairedCoordinateIndex) {
			std::cout << "(x";
			printNumberSign(-abscissaCoordinates[pairedCoordinateIndex]);
			std::cout << fabs(abscissaCoordinates[pairedCoordinateIndex]);
			std::cout << ")";
		}
	}
	std::cout << std::endl;
}

long double function(long double x) {
	return sin(x) + x;
}

int main() {
	std::cout << "How many coordinates do you want to enter?" << std::endl;
	size_t coordinatesNumber;
	std::cin >> coordinatesNumber;
	std::vector<long double> abscissaCoordinates(coordinatesNumber);
	std::cout << "Enter the coordinates: " << std::endl;
	for (size_t coordinateIndex = 0; coordinateIndex < coordinatesNumber; ++coordinateIndex) {
		std::cin >> abscissaCoordinates[coordinateIndex];
	}
	std::vector<Point> points = generatePointsByFunctionAndAbscissaCoordinates(function, abscissaCoordinates);
	std::vector<long double> newtonPolynomialCoefficients = findNewtonPolynomialCoefficients(points);
	std::cout << "Newton polynomial: " << std::endl;
	printPolynomial(newtonPolynomialCoefficients, abscissaCoordinates);
	long double checkCoordinate;
	std::cout << "Enter coordinate X*: " << std::endl;
	std::cin >> checkCoordinate;
	std::cout << "Newton polynomial value: " << std::endl;
	long double newtonPolynomialValue = calculatePolynomialValue(newtonPolynomialCoefficients, abscissaCoordinates, checkCoordinate);
	std::cout << newtonPolynomialValue << std::endl;
	std::cout << "Real function value: " << std::endl;
	long double realFunctionValue = function(checkCoordinate);
	std::cout << realFunctionValue << std::endl;
	std::cout << "Difference: " << std::endl;
	std::cout << fabs(newtonPolynomialValue - realFunctionValue) << std::endl;
}
