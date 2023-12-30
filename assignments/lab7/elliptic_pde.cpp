#include <iostream> 
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <tuple>
#include "../../matrix.h"
#include "../../matrix_operations.h"

#define EPS 0.000000000000001L

long double leftBoundaryTimeDerivativeFunction(long double y) {
  return 0;
}

long double rightBoundaryFunction(long double y) {
  return 1 - std::pow(y, 2);
}

long double bottomBoundaryTimeDerivativeFunction(long double x) {
  return 0;
}

long double topBoundaryFunction(long double x) {
  return std::pow(x, 2) - 1;
}

long double analyticalSolution(long double x, long double y) {
  return std::pow(x, 2) - std::pow(y, 2);
}

struct ProblemVariables {
  long long N1;
  long long N2;
  long double leftBoundary = 0;
  long double rightBoundary = 1;
  long double bottomBoundary = 0;
  long double topBoundary = 1;
  long double getStepX() const {
    return rightBoundary / (long double)N1;
  }
  long double getStepY() const {
    return topBoundary / (long double)N2;
  }
};

void getProblemVariables(ProblemVariables& problemVariables) {
  std::cout << "Enter N1: ";
  std::cin >> problemVariables.N1;
  std::cout << "Enter N2: ";
  std::cin >> problemVariables.N2;
}

std::vector<std::vector<long double>> getAnalyticSolution(const std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  size_t N1 = problemVariables.N1;
  size_t N2 = problemVariables.N2;
  long double h1 = problemVariables.getStepX();
  long double h2 = problemVariables.getStepY();
  for (size_t i = 0; i <= N1; ++i) {
    for (size_t j = 0; j <= N2; ++j) {
      u[i][j] = analyticalSolution(i * h1, j * h2);
    }
  }
  return u;
}

void fillBorderValues(std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  auto& u = finiteDifferenceGrid;
  for (size_t i = 0; i <= problemVariables.N1; ++i) {
    u[i][0] = u[i][1];
    u[i][problemVariables.N2] = topBoundaryFunction(i * problemVariables.getStepX());
  }
  for (size_t j = 0; j <= problemVariables.N2; ++j) {
    u[0][j] = u[1][j];
    u[problemVariables.N1][j] = rightBoundaryFunction(j * problemVariables.getStepY());
  }
} 

void fillInterpolation(std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  size_t N1 = problemVariables.N1;
  size_t N2 = problemVariables.N2;
  long double h1 = problemVariables.getStepX();
  long double h2 = problemVariables.getStepY();
  long double l1 = problemVariables.rightBoundary;
  long double l2 = problemVariables.topBoundary;
  auto& u = finiteDifferenceGrid;
  for (size_t i = 1; i < N1; ++i) {
    for (size_t j = 1; j < N2; ++j) {
      long double x = h1 * i;
      long double y = h2 * j;
      u[i][j] = rightBoundaryFunction(y) * (l2 - y) / ((l1 - x) + (l2 - y)) + topBoundaryFunction(x) * (l1 - x) / ((l1 - x) + (l2 - y));
    }
  }
}

long double getRootMeanSquare(const std::vector<std::vector<long double>>& scheme, const std::vector<std::vector<long double>>& exactSolution) {
  long double sum = 0;
  long double elementQuantity = 0;
  for (size_t i = 0; i < scheme.size(); ++i) { for (size_t j = 0; j < scheme[0].size(); ++j) {
      sum += std::pow(scheme[i][j] - exactSolution[i][j], 2);
      ++elementQuantity;
    }
  }
  return std::sqrt(sum / elementQuantity);
}

void printScheme(const std::vector<std::vector<long double>>& scheme) {
  for (const auto& line: scheme) {
    for (const auto& el: line) {
      std::cout << el << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

std::tuple<size_t, std::vector<std::vector<long double>>> solveLiebmannMethod(const std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  auto uPrev = finiteDifferenceGrid;
  long double h1 = problemVariables.getStepX();
  long double h2 = problemVariables.getStepY();
  fillInterpolation(u, problemVariables);
  fillBorderValues(u, problemVariables); 
  size_t iterationsCnt = 0;
  while (!(getRootMeanSquare(u, uPrev) <= EPS)) {
    ++iterationsCnt;
    uPrev = u;
    for (size_t i = 1; i < problemVariables.N1; ++i) {
      for (size_t j = 1; j < problemVariables.N2; ++j) {
        long double rightPart = 0;
        long double divider = 2 * (std::pow(h1, 2) + std::pow(h2, 2));
        u[i][j] = std::pow(h2, 2) / divider * (uPrev[i + 1][j] + uPrev[i - 1][j]) + std::pow(h1, 2) / divider * (uPrev[i][j - 1] + uPrev[i][j + 1]);
      }
    }
    fillBorderValues(u, problemVariables);
  }
  return make_tuple(iterationsCnt, u);
}

std::tuple<size_t, std::vector<std::vector<long double>>> solveSeidelMethod(const std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  auto u = finiteDifferenceGrid;
  auto uPrev = finiteDifferenceGrid;
  long double h1 = problemVariables.getStepX();
  long double h2 = problemVariables.getStepY();
  fillInterpolation(u, problemVariables);
  fillBorderValues(u, problemVariables); 
  size_t iterationsCnt = 0;
  while (!(getRootMeanSquare(u, uPrev) <= EPS)) {
    ++iterationsCnt;
    auto tmpPrev = uPrev;
    uPrev = u;
    for (size_t i = 1; i < problemVariables.N1; ++i) {
      for (size_t j = 1; j < problemVariables.N2; ++j) {
        long double rightPart = 0;
        long double divider = 2 * (std::pow(h1, 2) + std::pow(h2, 2));
        u[i][j] = std::pow(h2, 2) / divider * (u[i + 1][j] + u[i - 1][j]) + std::pow(h1, 2) / divider * (u[i][j - 1] + u[i][j + 1]);
      }
    }
    fillBorderValues(u, problemVariables);
  }
  return std::make_tuple(iterationsCnt, u);
}

std::tuple<size_t, std::vector<std::vector<long double>>> solveSimpleIterationsRelaxationsMethod(const std::vector<std::vector<long double>>& finiteDifferenceGrid, const ProblemVariables& problemVariables) {
  long double omega = 1.3;
  auto u = finiteDifferenceGrid;
  auto uPrev = finiteDifferenceGrid;
  long double h1 = problemVariables.getStepX();
  long double h2 = problemVariables.getStepY();
  fillInterpolation(u, problemVariables);
  fillBorderValues(u, problemVariables); 
  size_t iterationsCnt = 0;
  while (!(getRootMeanSquare(u, uPrev) <= EPS)) {
    ++iterationsCnt;
    auto tmpPrev = uPrev;
    uPrev = u;
    for (size_t i = 1; i < problemVariables.N1; ++i) {
      for (size_t j = 1; j < problemVariables.N2; ++j) {
        long double rightPart = 0;
        long double divider = 2 * (std::pow(h1, 2) + std::pow(h2, 2));
        u[i][j] = (1 - omega) * tmpPrev[i][j] + omega * (std::pow(h2, 2) / divider * (u[i + 1][j] + u[i - 1][j]) + std::pow(h1, 2) / divider * (u[i][j - 1] + u[i][j + 1]));
      }
    }
    fillBorderValues(u, problemVariables);
  }
  return std::make_tuple(iterationsCnt, u);
}

void outputToFileColumn(const std::string& fileName, size_t column, const ProblemVariables& problemVariables, 
                        const std::vector<std::vector<long double>>& scheme, const std::vector<std::vector<long double>>& solution) {
  std::ofstream timeSlice;
  timeSlice.open(fileName);
  for (size_t row = 0; row < scheme.size(); ++row) {
    timeSlice << problemVariables.getStepY() * row << " " << scheme[row][column] << std::endl;
  }
  timeSlice << std::endl;
  for (size_t row = 0; row < scheme.size(); ++row) {
    timeSlice << problemVariables.getStepY() * row << " " << solution[row][column] << std::endl;
  }
  timeSlice.close();
}

int main(int argc, char** argv) {
  ProblemVariables problemVariables;
  getProblemVariables(problemVariables);
  std::vector<std::vector<long double>> finiteDifferenceGrid(problemVariables.N1 + 1, std::vector<long double>(problemVariables.N2 + 1, 0));
  auto analytic = getAnalyticSolution(finiteDifferenceGrid, problemVariables);
  std::vector<std::vector<long double>> liebmannMethodResult;
  size_t liebmannIterationsQuantity;
  std::tie(liebmannIterationsQuantity, liebmannMethodResult) = solveLiebmannMethod(finiteDifferenceGrid, problemVariables);
  std::cout << "RMS Liebmann: " << getRootMeanSquare(liebmannMethodResult, analytic) << std::endl;
  std::cout << "Iterations Liebmann: " << liebmannIterationsQuantity << std::endl;
  outputToFileColumn("liebmann_data.gr", problemVariables.N2 / 2, problemVariables, liebmannMethodResult, analytic);
  std::vector<std::vector<long double>> seidelMethodResult; 
  size_t seidelIterationsQuantity;
  std::tie(seidelIterationsQuantity, seidelMethodResult)  = solveSeidelMethod(finiteDifferenceGrid, problemVariables);
  std::cout << "RMS Seidel: " << getRootMeanSquare(seidelMethodResult, analytic) << std::endl;
  std::cout << "Iterations Seidel: " << seidelIterationsQuantity << std::endl;
  outputToFileColumn("seidel_data.gr", problemVariables.N2 / 2, problemVariables, seidelMethodResult, analytic);
  std::vector<std::vector<long double>> simpleIterationsRelaxationsMethodResult;
  size_t simpleIterationsRelaxationsQuantity; 
  std::tie(simpleIterationsRelaxationsQuantity, simpleIterationsRelaxationsMethodResult) = solveSimpleIterationsRelaxationsMethod(finiteDifferenceGrid, problemVariables);
  std::cout << "RMS Simple Iterations Relaxations: " << getRootMeanSquare(seidelMethodResult, analytic) << std::endl;
  std::cout << "Iterations Simple Iterations Relaxations: " << simpleIterationsRelaxationsQuantity << std::endl;
  outputToFileColumn("relexations_data.gr", problemVariables.N2 / 2, problemVariables, simpleIterationsRelaxationsMethodResult, analytic);
  return 0;
}
