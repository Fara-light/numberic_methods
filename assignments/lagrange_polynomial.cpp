#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

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

std::vector<long double> findLagrangePolynomialCoefficients(const std::vector<Point>& points) {
	std::vector<long double> lagrangePolynomialCoefficients(points.size());
	for (size_t pointIndex = 0; pointIndex < points.size(); ++pointIndex) {
		Point currentPoint = points[pointIndex];
		long double pointW = 1;
		for (size_t pairedPointIndex = 0; pairedPointIndex < points.size(); ++pairedPointIndex) {
			if (pointIndex == pairedPointIndex) {
				continue;
			}
			Point currentPairedPoint = points[pairedPointIndex];
			pointW *= (currentPoint.x - currentPairedPoint.x);
		}
		long double pointCoefficient = currentPoint.y / pointW;
		lagrangePolynomialCoefficients[pointIndex] = pointCoefficient;
	}
	return lagrangePolynomialCoefficients;
}

long double calculatePolynomialValue(const std::vector<long double>& coefficients, const std::vector<long double>& abscissaCoordinates, long double x) {
	if (coefficients.size() != abscissaCoordinates.size()) {
		throw std::invalid_argument("printLagrangePolynomial, wrong arguments");
	}
	long double polynomialValue = 0;
	for(size_t pointIndex = 0; pointIndex < abscissaCoordinates.size(); ++pointIndex) {
		long double currentMonoidValue = coefficients[pointIndex];
		for (size_t pairedPointIndex = 0; pairedPointIndex < abscissaCoordinates.size(); ++pairedPointIndex) {
			if (pointIndex == pairedPointIndex) {
				continue;
			}
			long double pairedPoint = abscissaCoordinates[pairedPointIndex];
			currentMonoidValue *= (x - pairedPoint);
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
		throw std::invalid_argument("printLagrangePolynomial, wrong arguments");
	}
	for (size_t pointIndex = 0; pointIndex < abscissaCoordinates.size(); ++pointIndex) {
		if (pointIndex != 0 && coefficients[pointIndex] < 0) {
			printNumberSign(coefficients[pointIndex]);
		} 
		std::cout << std::fabs(coefficients[pointIndex]);
		for (size_t pairedPointIndex = 0; pairedPointIndex < abscissaCoordinates.size(); ++pairedPointIndex) {
			if (pointIndex == pairedPointIndex) {
				continue;
			}
			std::cout << "(x";
			printNumberSign(-abscissaCoordinates[pairedPointIndex]);
			std::cout << std::fabs(abscissaCoordinates[pairedPointIndex]);
			std::cout << ")";
		}
	}
	std::cout << std::endl;
}

long double function(long double x) {
	return sin(x) + x;
}

void lagrangeSection() {
	std::cout << "How many coordinates do you want to enter?" << std::endl;
	size_t coordinatesNumber;
	std::cin >> coordinatesNumber;
	std::vector<long double> abscissaCoordinates(coordinatesNumber);
	std::cout << "Enter the coordinates: " << std::endl;
	for (size_t coordinateIndex = 0; coordinateIndex < coordinatesNumber; ++coordinateIndex) {
		std::cin >> abscissaCoordinates[coordinateIndex];
	}
	std::vector<Point> points = generatePointsByFunctionAndAbscissaCoordinates(function, abscissaCoordinates);
	std::vector<long double> lagrangePolynomialCoefficients = findLagrangePolynomialCoefficients(points);
	std::cout << "Lagrange polynomial: " << std::endl;
	printPolynomial(lagrangePolynomialCoefficients, abscissaCoordinates);
	long double checkCoordinate;
	std::cout << "Enter coordinate X*: " << std::endl;
	std::cin >> checkCoordinate;
	std::cout << "Lagrange polynomial value: " << std::endl;
	long double lagrangePolynomialValue = calculatePolynomialValue(lagrangePolynomialCoefficients, abscissaCoordinates, checkCoordinate);
	std::cout << lagrangePolynomialValue << std::endl;
	std::cout << "Real function value: " << std::endl;
	long double realFunctionValue = function(checkCoordinate);
	std::cout << realFunctionValue << std::endl;
	std::cout << "Difference: " << std::endl;
	std::cout << fabs(lagrangePolynomialValue - realFunctionValue) << std::endl;
}

int main() {
	lagrangeSection();
	return 0;
}
