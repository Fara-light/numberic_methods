#include <iostream> 
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <tuple>
#include "../../matrix.h"
#include "../../matrix_operations.h"

struct ProblemVariables {
  long long I;
  long long J;
  long long T;
  long double mu1; 
  long double mu2;
  long double a;
  long double timeBoundary;
  long double getBoundaryX() const {
    return M_PI / 2 * mu1;
  }
  long double getBoundaryY() const {
    return M_PI / 2 * mu2;
  }
  long double getStepX() const {
    return this->getBoundaryX() / (long double)this->I;
  }
  long double getStepY() const {
    return this->getBoundaryY() / (long double)this->J;
  }
  long double getStepT() const {
    return timeBoundary / (long double)this->T;
  }
};

long double xLeftBoundaryFunction(long double x, long double t, const ProblemVariables& pv) {
  return std::cos(pv.mu1 * x) * std::exp(-(std::pow(pv.mu1, 2) + std::pow(pv.mu2, 2)) * pv.a * t);
}

long double xRightBoundaryFunction(long double x, long double t, const ProblemVariables& pv) {
  return 0;
}

long double yLeftBoundaryFunction(long double y, long double t, const ProblemVariables& pv) {
  return std::cos(pv.mu2 * y) * std::exp(-(std::pow(pv.mu1, 2) + std::pow(pv.mu2, 2)) * pv.a * t);
}

long double yRightBoundaryFunction(long double y, long double t, const ProblemVariables& pv) {
  return 0;
}

long double timeBoundary(long double x, long double y, const ProblemVariables& pv) {
  return std::cos(pv.mu1 * x) * std::cos(pv.mu2 * y);
}

long double analyticSolutionFunction(long double x, long double y, long double t, const ProblemVariables& pv) {
  return std::cos(pv.mu1 * x) * std::cos(pv.mu2 * y) * std::exp(-(std::pow(pv.mu1, 2) + std::pow(pv.mu2, 2)) * pv.a * t);
}

void getProblemVariables(ProblemVariables& pv) {
  std::cout << "Enter I: ";
  std::cin >> pv.I;
  std::cout << "Enter J: ";
  std::cin >> pv.J;
  std::cout << "Enter T: ";
  std::cin >> pv.T;
  std::cout << "Enter time boudnary: ";
  std::cin >> pv.timeBoundary;
  std::cout << "Enter a: ";
  std::cin >> pv.a;
  std::cout << "Enter mu1: ";
  std::cin >> pv.mu1;
  std::cout << "Enter mu2: ";
  std::cin >> pv.mu2;
}

void fillKnownValues(std::vector<std::vector<std::vector<long double>>>& u, const ProblemVariables& pv, bool halfTime = false) {
  for (size_t i = 0; i <= pv.I; ++i) {
    for (size_t k = 0; k < u[i][0].size(); ++k) {
      long double t = k * pv.getStepT();
      if (halfTime) {
        t /= 2;
      }
      u[i][0][k] = xLeftBoundaryFunction(i * pv.getStepX(), t, pv);
      u[i][pv.J][k] = xRightBoundaryFunction(i * pv.getStepX(), t, pv);
    }
  }
  for (size_t j = 0; j <= pv.J; ++j) {
    for (size_t k = 0; k < u[0][j].size(); ++k) {
      long double t = k * pv.getStepT();
      if (halfTime) {
        t /= 2;
      }
      u[0][j][k] = yLeftBoundaryFunction(j * pv.getStepY(), t, pv);
      u[pv.I][j][k] = yRightBoundaryFunction(j * pv.getStepY(), t, pv);
    }
  }
  for (size_t i = 0; i <= pv.I; ++i) {
    for (size_t j = 0; j <= pv.J; ++j) {
      u[i][j][0] = timeBoundary(i * pv.getStepX(), j * pv.getStepY(), pv);
    }
  }
}

#define DELETE_THIS 5

std::vector<std::vector<std::vector<long double>>> getAnalyticSolution(const ProblemVariables& pv) {
  std::vector<std::vector<std::vector<long double>>> u(pv.I + 1, std::vector<std::vector<long double>>(pv.J + 1, std::vector<long double>(pv.T + 1, 0)));
  for (size_t k = 0; k <= pv.T; ++k) {
    for (size_t i = 0; i <= pv.I; ++i) {
      for (size_t j = 0; j <= pv.J; ++j) {
        u[i][j][k] = analyticSolutionFunction(i * pv.getStepX(), j * pv.getStepY(), k * pv.getStepT(), pv);
      }
    }
  }
  return u;
}

std::vector<std::vector<std::vector<long double>>> solveAlternatingDirectionsMethod(const ProblemVariables& pv) {
  std::vector<std::vector<std::vector<long double>>> v(pv.I + 1, std::vector<std::vector<long double>>(pv.J + 1, std::vector<long double>(pv.T * 2 + 1 + DELETE_THIS, 0)));
  std::vector<std::vector<std::vector<long double>>> u(pv.I + 1, std::vector<std::vector<long double>>(pv.J + 1, std::vector<long double>(pv.T + 1, 0)));
  fillKnownValues(v, pv, true);
  long double h1 = pv.getStepX();
  long double h2 = pv.getStepY();
  long double tau = pv.getStepT();
  long double a = pv.a;
  for (size_t k = 0; k < pv.T; ++k) {
    for (size_t j = 1; j < pv.J; ++j) {
      Matrix<long double> A(pv.I - 1, 3);
      Matrix<long double> B(pv.I - 1, 1);
      for (size_t i = 1; i < pv.I; ++i) {
        if (i == 1) {
          A[i - 1][0] = 1;
          A[i - 1][1] = -(2 + 2 * std::pow(h1, 2) / (a * tau));
          A[i - 1][2] = 0;
          B[i - 1][0] = -std::pow(h1, 2) / std::pow(h2, 2) * (v[i][j + 1][k * 2] - 2 * v[i][j][k * 2] + v[i][j - 1][k * 2]) - 2 * std::pow(h1, 2) / (a * tau) * v[i][j][k * 2] - v[i - 1][j][k * 2 + 1];
        } else if (i + 1 == pv.I) {
          A[i - 1][0] = 0;
          A[i - 1][1] = -(2 + 2 * std::pow(h1, 2) / (a * tau));
          A[i - 1][2] = 1;
          B[i - 1][0] = -std::pow(h1, 2) / std::pow(h2, 2) * (v[i][j + 1][k * 2] - 2 * v[i][j][k * 2] + v[i][j - 1][k * 2]) - 2 * std::pow(h1, 2) / (a * tau) * v[i][j][k * 2] - v[i + 1][j][k * 2 + 1];
        } else {
          A[i - 1][0] = 1;
          A[i - 1][1] = -(2 + 2 * std::pow(h1, 2) / (a * tau));
          A[i - 1][2] = 1;
          B[i - 1][0] = -std::pow(h1, 2) / std::pow(h2, 2) * (v[i][j + 1][k * 2] - 2 * v[i][j][k * 2] + v[i][j - 1][k * 2]) - 2 * std::pow(h1, 2) / (a * tau) * v[i][j][k * 2];
        }
      }
      auto solution = solveTridiagonalMatrixSystemOfLineralEquations(A, B);
      for (size_t row = 0; row < pv.I - 1; ++row) {
        v[row + 1][j][k * 2 + 1] = solution[row][0];
      }
    }
    for (size_t i = 1; i < pv.I; ++i) {
      Matrix<long double> A(pv.J - 1, 3);
      Matrix<long double> B(pv.J - 1, 1);
      for (size_t j = 1; j < pv.J; ++j) {
        if (j == 1) {
          A[j - 1][0] = 1;
          A[j - 1][1] = -(2 + 2 * std::pow(h2, 2) / (a * tau));
          A[j - 1][2] = 0;
          B[j - 1][0] = -std::pow(h2, 2) / std::pow(h1, 2) * (v[i + 1][j][k * 2 + 1] - 2 * v[i][j][k * 2 + 1] + v[i - 1][j][k * 2 + 1]) - 2 * std::pow(h2, 2) / (a * tau) * v[i][j][k * 2 + 1] - v[i][j - 1][k * 2 + 2];
        } else if (j + 1 == pv.J) {
          A[j - 1][0] = 0;
          A[j - 1][1] = -(2 + 2 * std::pow(h2, 2) / (a * tau));
          A[j - 1][2] = 1;
          B[j - 1][0] = -std::pow(h2, 2) / std::pow(h1, 2) * (v[i + 1][j][k * 2 + 1] - 2 * v[i][j][k * 2 + 1] + v[i - 1][j][k * 2 + 1]) - 2 * std::pow(h2, 2) / (a * tau) * v[i][j][k * 2 + 1] - v[i][j + 1][k * 2 + 2];
        } else {
          A[j - 1][0] = 1;
          A[j - 1][1] = -(2 + 2 * std::pow(h2, 2) / (a * tau));
          A[j - 1][2] = 1;
          B[j - 1][0] = -std::pow(h2, 2) / std::pow(h1, 2) * (v[i + 1][j][k * 2 + 1] - 2 * v[i][j][k * 2 + 1] + v[i - 1][j][k * 2 + 1]) - 2 * std::pow(h2, 2) / (a * tau) * v[i][j][k * 2 + 1];
        }
      }
      auto solution = solveTridiagonalMatrixSystemOfLineralEquations(A, B);
      for (size_t row = 0; row < pv.J - 1; ++row) {
        v[i][row + 1][k * 2 + 2] = solution[row][0];
      }
    }
  }
  for (size_t k = 0; k <= pv.T; ++k) {
    for (size_t i = 0; i <= pv.I; ++i) {
      for (size_t j = 0; j <= pv.J; ++j) {
        u[i][j][k] = v[i][j][k * 2];
      }
    }
  }
  return u;
}

std::vector<std::vector<std::vector<long double>>> solvePartialStepsMethod(const ProblemVariables& pv) {
  std::vector<std::vector<std::vector<long double>>> v(pv.I + 1, std::vector<std::vector<long double>>(pv.J + 1, std::vector<long double>(pv.T * 2 + 1 + DELETE_THIS, 0)));
  std::vector<std::vector<std::vector<long double>>> u(pv.I + 1, std::vector<std::vector<long double>>(pv.J + 1, std::vector<long double>(pv.T + 1, 0)));
  fillKnownValues(v, pv, true);
  long double h1 = pv.getStepX();
  long double h2 = pv.getStepY();
  long double tau = pv.getStepT();
  long double a = pv.a;
  for (size_t k = 0; k < pv.T; ++k) {
    for (size_t j = 1; j < pv.J; ++j) {
      Matrix<long double> A(pv.I - 1, 3);
      Matrix<long double> B(pv.I - 1, 1);
      for (size_t i = 1; i < pv.I; ++i) {
        if (i == 1) {
          A[i - 1][0] = 0;
          A[i - 1][1] = -(2 + std::pow(h1, 2) / (a * tau));
          A[i - 1][2] = 1;
          B[i - 1][0] = -std::pow(h1, 2) / (a * tau) * v[i][j][k * 2] - v[i - 1][j][k * 2 + 1];
        } else if (i + 1 == pv.I) {
          A[i - 1][0] = 1;
          A[i - 1][1] = -(2 + std::pow(h1, 2) / (a * tau));
          A[i - 1][2] = 0;
          B[i - 1][0] = -std::pow(h1, 2) / (a * tau) * v[i][j][k * 2] - v[i + 1][j][k * 2 + 1];
        } else {
          A[i - 1][0] = 1;
          A[i - 1][1] = -(2 + std::pow(h1, 2) / (a * tau));
          A[i - 1][2] = 1;
          B[i - 1][0] = -std::pow(h1, 2) / (a * tau) * v[i][j][k * 2];
        }
      }
      auto solution = solveTridiagonalMatrixSystemOfLineralEquations(A, B);
      for (size_t row = 0; row < pv.I - 1; ++row) {
        v[row + 1][j][k * 2 + 1] = solution[row][0];
      }
    }
    for (size_t i = 1; i < pv.I; ++i) {
      Matrix<long double> A(pv.J - 1, 3);
      Matrix<long double> B(pv.J - 1, 1);
      for (size_t j = 1; j < pv.J; ++j) {
        if (j == 1) {
          A[j - 1][0] = 0;
          A[j - 1][1] = -(2 + std::pow(h2, 2) / (a * tau));
          A[j - 1][2] = 1;
          B[j - 1][0] = -std::pow(h2, 2) / (a * tau) * v[i][j][k * 2 + 1] - v[i][j - 1][k * 2 + 2];
        } else if (j + 1 == pv.J) {
          A[j - 1][0] = 1;
          A[j - 1][1] = -(2 + std::pow(h2, 2) / (a * tau));
          A[j - 1][2] = 0;
          B[j - 1][0] = -std::pow(h2, 2) / (a * tau) * v[i][j][k * 2 + 1] - v[i][ + 1][k * 2 + 2];
        } else {
          A[j - 1][0] = 1;
          A[j - 1][1] = -(2 + std::pow(h2, 2) / (a * tau));
          A[j - 1][2] = 1;
          B[j - 1][0] = -std::pow(h2, 2) / (a * tau) * v[i][j][k * 2 + 1];
        }
      }
      auto solution = solveTridiagonalMatrixSystemOfLineralEquations(A, B);
      for (size_t row = 0; row < pv.J - 1; ++row) {
        v[i][row + 1][k * 2 + 2] = solution[row][0];
      }
    }
  }
  for (size_t k = 0; k <= pv.T; ++k) {
    for (size_t i = 0; i <= pv.I; ++i) {
      for (size_t j = 0; j <= pv.J; ++j) {
        u[i][j][k] = v[i][j][k * 2];
      }
    }
  }
  return u;
}

long double getRootMeanSquare(const std::vector<std::vector<std::vector<long double>>>& scheme, const std::vector<std::vector<std::vector<long double>>>& exactSolution) {
  long double sum = 0;
  long double elementQuantity = 0;
  for (size_t i = 0; i < scheme.size(); ++i) {
    for (size_t j = 0; j < scheme[i].size(); ++j) {
      for (size_t k = 0; k < scheme[i][j].size(); ++k) {
        sum += std::pow(scheme[i][j][k] - exactSolution[i][j][k], 2);
        ++elementQuantity;
      }
    }
  }
  return std::sqrt(sum / elementQuantity);
}

void streamWriteTimeSlice(std::ofstream& os,
                          const ProblemVariables& pv,
                          std::vector<std::vector<std::vector<long double>>> scheme) {
  for (size_t i = 0; i <= pv.I; ++i) {
    for (size_t j = 0; j < pv.J; ++j) {
      os << i * pv.getStepX() << " " << j * pv.getStepY() << " " << scheme[i][j][pv.T] << std::endl;
    }
  }
}

void streamWriteCoordinateSlice(std::ofstream& os,
                          const ProblemVariables& pv,
                          std::vector<std::vector<std::vector<long double>>> scheme) {
  for (size_t j = 0; j <= pv.J; ++j) {
    for (size_t k = 0; k <= pv.T; ++k) {
      os << j * pv.getStepY() << " " << k * pv.getStepT() << " " << scheme[pv.I / 2][j][k] << std::endl;
    }
  }
}

void generateGnuplotFile(const std::vector<std::vector<std::vector<long double>>>& numeric, 
                         const std::vector<std::vector<std::vector<long double>>>& analytic,
                         const ProblemVariables& pv,
                         std::string baseFileName) {
  std::ofstream timeSlice;
  timeSlice.open(".time_slice_" + baseFileName + ".gr");
  streamWriteTimeSlice(timeSlice, pv, numeric);
  timeSlice << std::endl;
  streamWriteTimeSlice(timeSlice, pv, analytic);
  timeSlice.close();
  std::ofstream xSlice;
  xSlice.open(".x_slice_" + baseFileName + ".gr");
  streamWriteCoordinateSlice(xSlice, pv, numeric);
  xSlice << std::endl;
  streamWriteCoordinateSlice(xSlice, pv, analytic);
  xSlice.close();
}

int main(int argc, char** argv) {
  ProblemVariables pv;
  getProblemVariables(pv);
  auto partialMethodSolution = solvePartialStepsMethod(pv);
  auto alternatingMethodSolution = solveAlternatingDirectionsMethod(pv);
  auto analyticSolution = getAnalyticSolution(pv);
  std::cout << "Partial method solution RMS: " << getRootMeanSquare(partialMethodSolution, analyticSolution) << std::endl;
  std::cout << "Alternating method solution RMS: " << getRootMeanSquare(alternatingMethodSolution, analyticSolution) << std::endl;
  generateGnuplotFile(partialMethodSolution, analyticSolution, pv, "partial_method");
  generateGnuplotFile(alternatingMethodSolution, analyticSolution, pv, "alternating_method");
  return 0;
}
