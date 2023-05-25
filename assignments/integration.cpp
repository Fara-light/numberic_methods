#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <iomanip>

long double integrateRectangleMethod(long double (*function)(long double), long double left, long double right, long double h) {
	long double result = 0;
	while (left < right) {
		left += h;
		result += h * function((2 * left - h) / 2);
	}
	return result;
}

long double integrateTrapezoidMethod(long double (*function)(long double), long double left, long double right, long double h) {
	long double result = 0;
	while (left < right) {
		result += (function(left) + function(left + h)) * h / 2;
		left += h;
	}
	return result;
}

long double integrateSimpsonMethod(long double (*function)(long double), long double left, long double right, long double h) {
	long double previousX = left;
	long double currentX = left + h;
	long double result = h / 6 * (function(previousX) + 4 * function((previousX + currentX) / 2) + function(currentX));
	while (currentX < right) {
		previousX = currentX;
		currentX += h;
		result += h / 6 * (function(previousX) + 4 * function((previousX + currentX) / 2) + function(currentX));
	}
	return result;
}

// RungeRonbergRichardson

long double RungeRonbergRichardsonMethod(long double integralH1, long double integralH2, long double p) {
	return integralH1 + (integralH1 - integralH2) / (std::pow(2, p) - 1);
}

long double integrateRectangleMethodRRR(long double (*function)(long double), long double left, long double right, long double h) {
	long double integralH1 = integrateRectangleMethod(function, left, right, h);
	long double integralH2 = integrateRectangleMethod(function, left, right, h / 2);
	return RungeRonbergRichardsonMethod(integralH2, integralH1, 2);
}

long double integrateTrapezoidMethodRRR(long double (*function)(long double), long double left, long double right, long double h) {
	long double integralH1 = integrateTrapezoidMethod(function, left, right, h);
	long double integralH2 = integrateTrapezoidMethod(function, left, right, h / 2);
	return RungeRonbergRichardsonMethod(integralH2, integralH1, 2);
}

long double integrateSimpsonMethodRRR(long double (*function)(long double), long double left, long double right, long double h) {
	long double integralH1 = integrateSimpsonMethod(function, left, right, h);
	long double integralH2 = integrateSimpsonMethod(function, left, right, h / 2);
	return RungeRonbergRichardsonMethod(integralH2, integralH1, 4);
}

long double function(long double x) {
	return x / (pow(x, 3) + 8);
}

int main() {
	// stupid assigment makes me hard-code all the data
	long double left, right, h = 0.5, h2 = 0.25;
	std::cout << "Enter left boundary: " << std::endl;
	std::cin >> left;
	std::cout << "Enter right boundary: " << std::endl;
	std::cin >> right;
	std::cout << "h1 = 0.5, h2 = 0.25" << std::endl;
	std::cout << std::setprecision(15);
	long double rectangleMethod = integrateRectangleMethod(function, left, right, h);
	long double rectangleMethod2 = integrateRectangleMethod(function, left, right, h2);
	long double rectangleMethodRRR = integrateRectangleMethodRRR(function, left, right, h);
	std::cout << "RectangleMethod, h = 0.5: " << rectangleMethod << std::endl;
	std::cout << "RectangleMethod, h = 0.25: " << rectangleMethod2 << std::endl;
	std::cout << "RectangleMethod RRR: " << rectangleMethodRRR << std::endl;
	std::cout << "Difference: " << rectangleMethod - rectangleMethodRRR << std::endl;
	long double trapezoidMethod = integrateTrapezoidMethod(function, left, right, h);
	long double trapezoidMethod2 = integrateTrapezoidMethod(function, left, right, h2);
	long double trapezoidMethodRRR = integrateTrapezoidMethodRRR(function, left, right, h);
	std::cout << "Trapezoid method, h = 0.5: " << trapezoidMethod << std::endl;
	std::cout << "Trapezoid method, h = 0.25: " << trapezoidMethod2 << std::endl;
	std::cout << "Trapezoid method RRR: " << trapezoidMethodRRR << std::endl;
	std::cout << "Difference: " << trapezoidMethod - trapezoidMethodRRR << std::endl;
	long double simpsonMethod = integrateSimpsonMethod(function, left, right, h);
	long double simpsonMethod2 = integrateSimpsonMethod(function, left, right, h2);
	long double simpsonMethodRRR = integrateSimpsonMethodRRR(function, left, right, h);
	std::cout << "Simpson method, h = 0.5: " << simpsonMethod << std::endl;
	std::cout << "Simpson method, h = 0.25: " << simpsonMethod2 << std::endl;
	std::cout << "Simpson method RRR: " << simpsonMethodRRR << std::endl;
	std::cout << "Difference: " << simpsonMethod - simpsonMethodRRR << std::endl;
}
