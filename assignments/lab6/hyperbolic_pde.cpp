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
  second = 2,
  third = 3
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
    } else if (order == ApproximationOrder::third) {
      u[k + 1][0] = (u[k + 1][0] + 1 / sigma * u[k][0] - 1 / (2 * sigma) * u[k - 1][0]) / (1 + problemVariables.getSpatialStep() + 1 / (2 * sigma));
      u[k + 1][problemVariables.N] = (u[k + 1][problemVariables.N - 1] + 1 / sigma * u[k][problemVariables.N] - 1 / (2 * sigma) * u[k - 1][problemVariables.N]) / (1 - problemVariables.getSpatialStep() + 1 / (2 * sigma));
    }
  }
  return u;
}

std::vector<std::vector<long double>> getImplicitFiniteDifferenceScheme(const std::vector<std::vector<long double>> finiteDifferenceGrid, const ProblemVariables& problemVariables, const ApproximationOrder order) {
  auto u = finiteDifferenceGrid;
  long double omega = -std::pow(problemVariables.a, 2) * std::pow(problemVariables.getTimeStep(), 2);
  long double sigma = std::pow(problemVariables.a, 2) * std::pow(problemVariables.getTimeStep(), 2) / std::pow(problemVariables.getSpatialStep(), 2);
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
        } else if (order == ApproximationOrder::second) {
          a = 0;
          b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega + 4 * (-omega) / (3 + 2 * problemVariables.getSpatialStep());
          c = -omega - (-omega) / (3 + 2 * problemVariables.getSpatialStep());
        } else if (order == ApproximationOrder::third) {
          a = 0;
          b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega + (-omega) / (1 + problemVariables.getSpatialStep() + 1 / (2 * sigma));
          c = -omega;
          d = d - (-omega) * (1 / sigma * u[k][0] - 1 / (2 * sigma) * u[k - 1][0]) / (1 + problemVariables.getSpatialStep() + 1 / (2 * sigma));
        }
      } if (j < problemVariables.N - 1 && j > 1) {
        a = -omega;
        b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega;
        c = -omega;
      } else if (j == problemVariables.N - 1) {
        c = 0;
        if (order == ApproximationOrder::first) {
          a = -omega;
          b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega + -omega / (1 - problemVariables.getSpatialStep());
          c = 0;
        } else if (order == ApproximationOrder::second) {
          a = -omega - (-omega) / (3 - problemVariables.getSpatialStep());
          b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega + 4 * (-omega) / (3 - problemVariables.getSpatialStep()); 
          c = 0;
        } else if (order == ApproximationOrder::third) {
          a = -omega;
          b = std::pow(problemVariables.getSpatialStep(), 2) + 2 * omega + (-omega) / (1 - problemVariables.getSpatialStep() + 1 / (2 * sigma));
          c = 0;
          d = d - (-omega) * (1 / sigma * u[k][problemVariables.N] - 1 / (2 * sigma) * u[k - 1][problemVariables.N])/(1 - problemVariables.getSpatialStep() + 1 / (2 * sigma));
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
    } else if (order == ApproximationOrder::third) {
      u[k + 1][0] = (u[k + 1][0] + 1 / sigma * u[k][0] - 1 / (2 * sigma) * u[k - 1][0]) / (1 + problemVariables.getSpatialStep() + 1 / (2 * sigma));
      u[k + 1][problemVariables.N] = (u[k + 1][problemVariables.N - 1] + 1 / sigma * u[k][problemVariables.N] - 1 / (2 * sigma) * u[k - 1][problemVariables.N]) / (1 - problemVariables.getSpatialStep() + 1 / (2 * sigma));
    }
  }
  return u;
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

void printScheme(const std::vector<std::vector<long double>>& grid) {
  for (const auto& line: grid) {
    for (const auto& el: line) {
      std::cout << el << " ";
    }
    std::cout << std::endl;
  }
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

int main() {
  ProblemVariables problemVariables;
  getProblemVariables(problemVariables);
  long double sigma = std::pow(problemVariables.a, 2) * std::pow(problemVariables.getTimeStep(), 2) / std::pow(problemVariables.getSpatialStep(), 2);
  std::vector<std::vector<long double>> finiteDifferenceGrid1(problemVariables.K + 1, std::vector<long double>(problemVariables.N + 1, 0));
  fillKnownGridValues(finiteDifferenceGrid1, problemVariables);
  fillFirstStepValues(finiteDifferenceGrid1, problemVariables, &firstStateFirstOrderApproximation);
  std::vector<std::vector<long double>> finiteDifferenceGrid2(problemVariables.K + 1, std::vector<long double>(problemVariables.N + 1, 0));
  fillKnownGridValues(finiteDifferenceGrid2, problemVariables);
  fillFirstStepValues(finiteDifferenceGrid2, problemVariables, &firstStateSecondOrderApproximation);
  auto analyticSolution = getAnalyticSolution(finiteDifferenceGrid1, problemVariables);
  std::vector<std::vector<long double>> explicitFiniteDifferenceScheme1_1, explicitFiniteDifferenceScheme1_2, explicitFiniteDifferenceScheme1_3;
  std::vector<std::vector<long double>> explicitFiniteDifferenceScheme2_1, explicitFiniteDifferenceScheme2_2, explicitFiniteDifferenceScheme2_3;
  if (sigma < 1) {
    explicitFiniteDifferenceScheme1_1 = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid1, problemVariables, ApproximationOrder::first); 
    std::cout << "Explicit, 1st order initial, 1st border. RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme1_1, analyticSolution) << std::endl;
    explicitFiniteDifferenceScheme1_2 = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid1, problemVariables, ApproximationOrder::second);
    std::cout << "Explicit, 1st order initial, 2d border. RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme1_2, analyticSolution) << std::endl;
    explicitFiniteDifferenceScheme1_3 = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid1, problemVariables, ApproximationOrder::third);
    std::cout << "Explicit, 1d order initial, 2d border - 2 points. RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme1_3, analyticSolution) << std::endl;
    explicitFiniteDifferenceScheme2_1 = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid2, problemVariables, ApproximationOrder::first); 
    std::cout << "Explicit, 2d order initial, 1st border. RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme2_1, analyticSolution) << std::endl;
    explicitFiniteDifferenceScheme2_2 = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid2, problemVariables, ApproximationOrder::second);
    std::cout << "Explicit, 2d order initial, 2d border. RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme2_2, analyticSolution) << std::endl;
    explicitFiniteDifferenceScheme2_3 = getExplicitFiniteDifferenceScheme(finiteDifferenceGrid2, problemVariables, ApproximationOrder::third);
    std::cout << "Explicit, 2d order initial, 2d border - 2 points. RMS: " << getRootMeanSquare(explicitFiniteDifferenceScheme2_3, analyticSolution) << std::endl;
  } else {
    std::cout << "Impossible to build explicit finite difference scheme: sigma >= 1" << std::endl;
  }
  auto implicitFiniteDifferenceScheme1_1 = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid1, problemVariables, ApproximationOrder::first); 
  std::cout << "Implicit, 1st order initial, 1st border. RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme1_1, analyticSolution) << std::endl;
  auto implicitFiniteDifferenceScheme1_2 = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid1, problemVariables, ApproximationOrder::second);
  std::cout << "Implicit, 1st order initial, 2d border. RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme1_2, analyticSolution) << std::endl;
  auto implicitFiniteDifferenceScheme1_3 = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid1, problemVariables, ApproximationOrder::third);
  std::cout << "Implicit, 1st order initial, 2d border - 2 points. RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme1_3, analyticSolution) << std::endl;
  auto implicitFiniteDifferenceScheme2_1 = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid2, problemVariables, ApproximationOrder::first);
  std::cout << "Implicit, 2d order initial, 1st border. RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme2_1, analyticSolution) << std::endl;
  auto implicitFiniteDifferenceScheme2_2 = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid2, problemVariables, ApproximationOrder::second);
  std::cout << "Implicit, 2d order initial, 2d border. RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme2_2, analyticSolution) << std::endl;
  auto implicitFiniteDifferenceScheme2_3 = getImplicitFiniteDifferenceScheme(finiteDifferenceGrid2, problemVariables, ApproximationOrder::third);
  std::cout << "Implicit, 2d order initial, 2d border - 2 points. RMS: " << getRootMeanSquare(implicitFiniteDifferenceScheme2_3, analyticSolution) << std::endl;
  outputToFileColumn(".result", problemVariables.N, problemVariables, explicitFiniteDifferenceScheme2_3, analyticSolution);
  return 0;
}
