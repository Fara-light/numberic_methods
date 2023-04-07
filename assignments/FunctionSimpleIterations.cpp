#include <iostream>
#include <cmath>
using namespace std;

long double phi(long double x) {
	return sqrt((pow(3, x) + 1) / 5);
}

long double dphi_dx(long double x) {
	return (sqrt(5 * pow(3, x) + 5) * pow(3, x) * log(3)) / (10 * pow(3, x) + 10);
}

std::pair<size_t, long double> simpleIterations(long double leftEnd, long double rightEnd, long double EPS) {
	long double q = fabs(dphi_dx(rightEnd));
	if (q >= 1) {
		throw::range_error("can't use simple iterations when q >= 1");
	}
	long double x = (leftEnd + rightEnd) / 2;
	long double prevX;
	size_t iterationsCount = 0;
	do {
		prevX = x;
		x = phi(x);
		++iterationsCount;
	} while (q * fabs(x - prevX) / (1 - x) >= EPS);
	return make_pair(iterationsCount, x);
}

int main() {
	long double leftEnd, rightEnd;
	std::cout << "Enter left end: ";
	std::cin >> leftEnd;
	std::cout << "Enter right end: ";
	std::cin >> rightEnd;
	long double EPS;
	std::cout << "Enter EPS: ";
	std::cin >> EPS;
	std::pair<size_t, long double> solution = simpleIterations(leftEnd, rightEnd, EPS);
	std::cout << "Simple iterations method finished in " << solution.first << " iterations" << std::endl;
	std::cout << "Ans: " << solution.second << std::endl;
	return 0;
}
