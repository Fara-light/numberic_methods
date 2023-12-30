#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../../matrix.h"
#include "../../matrix_operations.h"

long double leftBoundaryFunction(long double a, long double t) {
  return std::exp(-a * t);
}

long double rightBoundaryFunction(long double a, long double t) {
  return -std::exp(-a * t);
}

long double initialStateFunction(long double x) {
  return std::cos(x);
}

long double analyticalSolution(long double a, long double x, long double t) {
  return std::exp(-a * t) * std::cos(x);
}

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

long double getTime(size_t k, const ProblemVariables& problemVariables) {
  return k * problemVariables.getTimeStep();
}

long double getCoordinate(size_t j, const ProblemVariables& problemVariables) {
  return j * problemVariables.getSpatialStep();
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
  for (size_t k = 0; k < problemVariables.K + 1; ++k) {
    long double time = getTime(k, problemVariables);
    finiteDifferenceGrid[k][0] = leftBoundaryFunction(problemVariables.a, time);
    finiteDifferenceGrid[k][problemVariables.N] = rightBoundaryFunction(problemVariables.a, time);
 }
}

std::vector<std::vector<long double>> getExplicitFiniteDifferenceScheme(const std::vector<std::vector<long double>>& finiteDifferenceGrid,
                                                                        const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  for (size_t k = 1; k < problemVariables.K + 1; ++k) {
    for (size_t j = 1; j < problemVariables.N; ++j) {
      long double gamma = problemVariables.a * problemVariables.getTimeStep();
      long double sigma = (u[k - 1][j + 1] - 2 * u[k - 1][j] + u[k - 1][j - 1]) / std::pow(problemVariables.getSpatialStep(), 2);
      u[k][j] = gamma * sigma + u[k - 1][j];
    }
  }
  return u;
}

std::vector<std::vector<long double>> getImplicitFiniteDifferenceScheme(const std::vector<std::vector<long double>>& finiteDifferenceGrid,
                                                                        const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  long double sigma = problemVariables.a * problemVariables.getTimeStep() / std::pow(problemVariables.getSpatialStep(), 2);
  for (size_t k = 0; k < problemVariables.K; ++k) {
    Matrix<long double> threeDiagonalMatrix(problemVariables.N - 1, 3); 
    Matrix<long double> rightPartMatrix(problemVariables.N - 1, 1);
    for (size_t row = 0; row < problemVariables.N - 1; ++row) {
      long double a, b, c, d;
      a = sigma;
      b = -(1 + 2 * sigma);
      c = sigma;
      if (row == 0) {
        d = -(u[k][1] + sigma * u[k + 1][0]);
      } else if (row == problemVariables.N - 2) {
        d = -(u[k][problemVariables.N - 1] + sigma * u[k + 1][problemVariables.N]);
      } else {
        d = -u[k][row + 1];
      }
      threeDiagonalMatrix[row][0] = a;
      threeDiagonalMatrix[row][1] = b;
      threeDiagonalMatrix[row][2] = c;
      rightPartMatrix[row][0] = d;
    }
    auto solution = solveTridiagonalMatrixSystemOfLineralEquations(threeDiagonalMatrix, rightPartMatrix);
    for (size_t j = 1; j <= problemVariables.N - 1; ++j) {
      u[k + 1][j] = solution[j - 1][0];
    }
  }
  return u;
}

std::vector<std::vector<long double>> getCrankNicolsonScheme(const std::vector<std::vector<long double>>& finiteDifferenceGrid,
                                                             const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  long double theta = 0.5;
  for (size_t k = 1; k <= problemVariables.K; ++k) {
    Matrix<long double> threeDiagonalMatrix(problemVariables.N - 1, 3);
    Matrix<long double> rightPartMatrix(problemVariables.N - 1, 1);
    for (size_t row = 0; row < problemVariables.N - 1; ++row) {
      long double a, b, c, d;
      long double sigma = problemVariables.a * problemVariables.getTimeStep() / std::pow(problemVariables.getSpatialStep(), 2);
      a = -theta * sigma;
      b = 2 * theta * sigma + 1;
      c = -theta * sigma;
      if (row == 0) {
        d = u[k - 1][row + 1] + (1 - theta) * sigma * (u[k - 1][row + 2] - 2 * u[k - 1][row + 1] + u[k - 1][row]) + theta * sigma * u[k][0];
      } else if (row == problemVariables.N - 2) {
        d = u[k - 1][row + 1] + (1 - theta) * sigma * (u[k - 1][row + 2] - 2 * u[k - 1][row + 1] + u[k - 1][row]) + theta * sigma * u[k][problemVariables.N];
      } else {
        d = u[k - 1][row + 1] + (1 - theta) * sigma * (u[k - 1][row + 2] - 2 * u[k - 1][row + 1] + u[k - 1][row]);
      }
      threeDiagonalMatrix[row][0] = a;
      threeDiagonalMatrix[row][1] = b;
      threeDiagonalMatrix[row][2] = c;
      rightPartMatrix[row][0] = d;
    }
    auto solution = solveTridiagonalMatrixSystemOfLineralEquations(threeDiagonalMatrix, rightPartMatrix);
    for (size_t j = 1; j <= problemVariables.N - 1; ++j) {
      u[k][j] = solution[j - 1][0];
    }
  }
  return u;
}

std::vector<std::vector<long double>> getAnswer(const std::vector<std::vector<long double>>& finiteDifferenceGrid,
                                                const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  for (size_t k = 0; k < problemVariables.K + 1; ++k) {
    for (size_t j = 0; j < problemVariables.N + 1; ++j) {
      long double time = getTime(k, problemVariables);
      long double x = getCoordinate(j, problemVariables);
      u[k][j] = analyticalSolution(problemVariables.a, x, time);
    }
  }
  return u;
}

void printScheme(const std::vector<std::vector<long double>>& scheme) { 
  for (const auto& row: scheme) {
    for (const auto& el: row) {
      std::cout << el << " ";
    }
    std::cout << std::endl;
  }
}

long double getRootMeanSquare(const std::vector<std::vector<long double>>& scheme, const std::vector<std::vector<long double>>& exactSolution) {
  long double sum = 0;
  long double elementQuantity = 0;
  for (size_t i = 0; i < scheme.size(); ++i) {
    for (size_t j = 0; j < scheme[0].size(); ++j) {
      sum += std::pow(scheme[i][j] - exactSolution[i][j], 2);
      ++elementQuantity;
    }
  }
  return std::sqrt(sum / elementQuantity);
}

void outputToFileColumn(const std::string& fileName, size_t column, const ProblemVariables& problemVariables, 
                        const std::vector<std::vector<long double>>& scheme, const std::vector<std::vector<long double>>& solution) {
  std::ofstream timeSlice;
  timeSlice.open(fileName);
  for (size_t row = 0; row < scheme.size(); ++row) {
    timeSlice << problemVariables.getSpatialStep() * row << " " << scheme[row][column] << std::endl;
  }
  timeSlice << std::endl;
  for (size_t row = 0; row < scheme.size(); ++row) {
    timeSlice << problemVariables.getSpatialStep() * row << " " << solution[row][column] << std::endl;
  }
  timeSlice.close();
}

void clearFile(const std::string& fileName) {
  std::ofstream timeSlice;
  timeSlice.open(fileName);
  timeSlice.close();
}

int main() {
  ProblemVariables problemVariables;
  getProblemVariables(problemVariables);
  std::vector<std::vector<long double>> finiteDifferenceGrid(problemVariables.K + 1, std::vector<long double>(problemVariables.N + 1, 0));
  fillKnownGridValues(finiteDifferenceGrid, problemVariables);
  long double eta = problemVariables.getTimeStep() / std::pow(problemVariables.getSpatialStep(), 2) * problemVariables.a;
  auto solution = getAnswer(finiteDifferenceGrid, problemVariables);
  if (eta >= 0.5) {
    std::cout << "System is not stable\n";
    std::cout << "Won't execute explicit finite difference scheme\n";
    clearFile(".explicitDifferenceScheme");
  } else { 
    std::cout << "Explicit finite difference scheme: " << std::endl;
    auto explicitFiniteDifferenceScheme = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid, problemVariables);
    outputToFileColumn(".explicitDifferenceScheme", problemVariables.N, problemVariables, explicitFiniteDifferenceScheme, solution);
    std::cout << "RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme, solution) << std::endl;
  } 
  std::cout << "Implicit finite difference scheme: " << std::endl;
  auto implicitFiniteDifferenceScheme = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid, problemVariables);
  outputToFileColumn(".implicitDifferenceScheme", problemVariables.N, problemVariables, implicitFiniteDifferenceScheme, solution);
  std::cout << "RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme, solution) << std::endl;
  std::cout << "Crank-Nicolson method: " << std::endl;
  auto crankNicolsonScheme = getCrankNicolsonScheme(finiteDifferenceGrid, problemVariables);
  outputToFileColumn(".crankNicolsonDifferenceScheme", problemVariables.N, problemVariables, crankNicolsonScheme, solution);
  std::cout << "RMS: " << getRootMeanSquare(crankNicolsonScheme, solution) << std::endl;
  return 0;
}
