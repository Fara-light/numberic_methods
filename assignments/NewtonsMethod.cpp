#include <iostream>
#include <cmath>
#include <stdexcept>

long double f(long double x) {
	return pow(3, x) - 5 * pow(x, 2) + 1;
}

long double df_dx(long double x) {
	return -10 * x + pow(3, x) * log(3);
}

std::pair<size_t, long double> newtonMethod(long double leftEnd, long double rightEnd, long double EPS) {
	if (leftEnd >= rightEnd) {
		throw std::invalid_argument("left end it bigger, than right end");
	}
	long double x = leftEnd + (rightEnd - leftEnd) / 2;
	long double prevX;
	size_t iterationsCount = 0;
	do {
		long double fValue = f(x);
		prevX = x;
		x = x - fValue / df_dx(x);
		++iterationsCount;

	} while (std::fabs(x - prevX) >= EPS);
	return std::make_pair(iterationsCount, x);
}

int main() {
	long double leftEnd, rightEnd;
	long double EPS;
	std::cout << "Enter left end: ";
	std::cin >> leftEnd;
	std::cout << "Enter right end: ";
	std::cin >> rightEnd;
	std::cout << "Enter EPS: ";
	std::cin >> EPS;
	std::pair<size_t, long double> solution = newtonMethod(leftEnd, rightEnd, EPS);
	std::cout << "Newton method finished in " << solution.first << " iterations" << std::endl;
	std::cout << "Ans: " << solution.second << std::endl;
	return 0;
}
