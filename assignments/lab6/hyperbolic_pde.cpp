#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../../matrix.h"
#include "../../matrix_operations.h"

struct ProblemVariables {
  long double a;
  long long N;
  long long K;
  long double T;
  long double leftBoundary = 0;
  long double rightBoundary = M_PI;
  long double getTimeStep() const {
    return T / K;
  }
  long double getSpatialStep() const {
    return rightBoundary / N;
  }
};

long double analyticSolution(long double a, long double x, long double t) {
  return std::sin(x - a * t) + std::cos(x + a * t);
}

long double initialStateFunction(long double x) {
  return std::sin(x) + std::cos(x);
}

long double initialStateTimeDerivativeFunction(long double a, long double x) {
  return - a * (std::sin(x) + std::cos(x));
}

long double firstStateFirstOrderApproximation(long double a, long double x, const ProblemVariables& problemVariables) {
  return initialStateFunction(x) + initialStateTimeDerivativeFunction(a, x) * problemVariables.getTimeStep();
}

long double firstStateSecondOrderApproximation(long double a, long double x, const ProblemVariables& problemVariables) {
  return initialStateFunction(x) + initialStateTimeDerivativeFunction(a, x) * problemVariables.getTimeStep() + a * a * (-std::sin(x) - std::cos(x)) * std::pow(problemVariables.getTimeStep(), 2) / 2;
}

void getProblemVariables(ProblemVariables& problemVariables) { 
  std::cout << "Enter a: ";
  std::cin >> problemVariables.a;
  std::cout << "Enter N: ";
  std::cin >> problemVariables.N;
  std::cout << "Enter K: ";
  std::cin >> problemVariables.K;
  std::cout << "Enter T: ";
  std::cin >> problemVariables.T;
  std::cout << "Result: " << std::endl;
}

void fillKnownGridValues(std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  for (size_t j = 0; j < problemVariables.N + 1; ++j) {
    finiteDifferenceGrid[0][j] = initialStateFunction(j * problemVariables.getSpatialStep());
  }
}

void fillFirstStepValues(std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables,
                         long double (*firstStepFunction)(long double, long double, const ProblemVariables&)) {
  for (size_t j = 0; j < problemVariables.N + 1; ++j) {
    finiteDifferenceGrid[1][j] = firstStepFunction(problemVariables.a, j * problemVariables.getSpatialStep(), problemVariables);
  }
}

enum class ApproximationOrder {
  first = 1,
  second = 2
};

std::vector<std::vector<long double>> getAnalyticSolution(const std::vector<std::vector<long double>> finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  for (size_t k = 1; k <= problemVariables.K; ++k) {
    for (size_t j = 0; j <= problemVariables.N; ++j) {
      u[k][j] = analyticSolution(problemVariables.a, j * problemVariables.getSpatialStep(), k * problemVariables.getTimeStep());
    }
  }
  return u;
}

std::vector<std::vector<long double>> getExplicitFiniteDifferenceScheme(const std::vector<std::vector<long double>> finiteDifferenceGrid, const ProblemVariables& problemVariables, const ApproximationOrder order) {
  auto u = finiteDifferenceGrid;
  long double sigma = std::pow(problemVariables.a, 2) * std::pow(problemVariables.getTimeStep(), 2) / std::pow(problemVariables.getSpatialStep(), 2);
  for (size_t k = 1; k < problemVariables.K; ++k) {
    for (size_t j = 1; j < problemVariables.N; ++j) {
      u[k + 1][j] = sigma * u[k][j + 1] + 2 * (1 - sigma) * u[k][j] + sigma * u[k][j - 1] - u[k - 1][j];
    }
    if (order == ApproximationOrder::first) {
      u[k + 1][0] = u[k + 1][1] / (1 + problemVariables.getSpatialStep());
      u[k + 1][problemVariables.N] = u[k + 1][problemVariables.N - 1] / (1 - problemVariables.getSpatialStep());
    } else if (order == ApproximationOrder::second) {
      u[k + 1][0] = (4 * u[k + 1][1] - u[k + 1][2]) / (3 + 2 * problemVariables.getSpatialStep());
      u[k + 1][problemVariables.N] = (4 * u[k + 1][problemVariables.N - 1] - u[k + 1][problemVariables.N - 2]) / (3 - 2 * problemVariables.getSpatialStep());
    }
  }
  return u;
}

std::vector<std::vector<long double>> getImplicitFiniteDifferenceScheme(const std::vector<std::vector<long double>> finiteDifferenceGrid, const ProblemVariables& problemVariables, const ApproximationOrder order) {
  auto u = finiteDifferenceGrid;
  long double omega = -std::pow(problemVariables.a, 2) * std::pow(problemVariables.getTimeStep(), 2);
  Matrix<long double> threeDiagonalMatrix(problemVariables.N - 1, 3);
  Matrix<long double> rightPartMatrix(problemVariables.N - 1, 1);
  for (size_t k = 1; k < problemVariables.K; ++k) {
    for (size_t j = 1; j < problemVariables.N; ++j) {
      long double a, b, c, d;
      d = 2 * std::pow(problemVariables.getSpatialStep(), 2) * u[k][j] - std::pow(problemVariables.getSpatialStep(), 2) * u[k - 1][j];
      if (j == 1) {
        a = 0;
        if (order == ApproximationOrder::first) {
          a = 0;
          b = -omega / (1 + problemVariables.getSpatialStep()) + std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega;
          c = -omega;
        }
      } if (j == problemVariables.N - 1) {
        a = -omega;
        b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega;
        c = -omega;
      } else {
        if (order == ApproximationOrder::second) {
          a = -omega;
          b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega + -omega / (1 - problemVariables.getSpatialStep());
          c = 0;
        }
      }
      threeDiagonalMatrix[j - 1][0] = a;
      threeDiagonalMatrix[j - 1][1] = b;
      threeDiagonalMatrix[j - 1][2] = c;
      rightPartMatrix[j - 1][0] = d;
    }
    auto solution = solveTridiagonalMatrixSystemOfLineralEquations(threeDiagonalMatrix, rightPartMatrix);
    for (size_t j = 1; j < problemVariables.N; ++j) {
      u[k + 1][j] = solution[j - 1][0];
    }
    if (order == ApproximationOrder::first) {
      u[k + 1][0] = u[k + 1][1] / (1 + problemVariables.getSpatialStep());
      u[k + 1][problemVariables.N] = u[k + 1][problemVariables.N - 1] / (1 - problemVariables.getSpatialStep());
    } else if (order == ApproximationOrder::second) {
      u[k + 1][0] = (4 * u[k + 1][1] - u[k + 1][2]) / (3 + 2 * problemVariables.getSpatialStep());
      u[k + 1][problemVariables.N] = (4 * u[k + 1][problemVariables.N - 1] - u[k + 1][problemVariables.N - 2]) / (3 - 2 * problemVariables.getSpatialStep());
    }
  }
  return u;
}

int main() {
  ProblemVariables problemVariables;
  getProblemVariables(problemVariables);
  std::vector<std::vector<long double>> finiteDifferenceGrid(problemVariables.K + 1, std::vector<long double>(problemVariables.N + 1, 0));
  fillKnownGridValues(finiteDifferenceGrid, problemVariables);
  fillFirstStepValues(finiteDifferenceGrid, problemVariables, &firstStateSecondOrderApproximation);
  auto explicitFiniteDifferenceScheme = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid, problemVariables, ApproximationOrder::first);
  auto analyticSolution = getAnalyticSolution(finiteDifferenceGrid, problemVariables);
  std::cout << "Explicit finite difference scheme, first order border approximation:\n";
  for (const auto line: explicitFiniteDifferenceScheme) {
    for (const auto p: line) {
      std::cout << p << " ";
    }
    std::cout << "\n";
  }
  std::cout << "Analytic solution:\n";
  for (const auto line: analyticSolution) {
    for (const auto p: line) {
      std::cout << p << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
