#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <assert.h>

template <typename Type> class MatrixView {
private:
  std::vector<Type> &v;
  MatrixView(const MatrixView &);
  MatrixView &operator=(const MatrixView &);

public:
  const size_t N, M;
  MatrixView(std::vector<Type> &v, size_t N, size_t M) : v(v), N(N), M(M) {
    assert(v.size() / N == M);
  }
  Type &set(size_t i, size_t j) { return v[i + N * j]; }
  const Type &get(size_t i, size_t j) { return v[i + N * j]; }
  Type &set(size_t n) { return v[n]; }
  const Type &get(size_t n) { return v[n]; }
};

double ParticularSolution(double x, double y) {
  return sin(2 * M_PI * x) * sinh(2 * M_PI * y);
}

double NormL2(const std::vector<double> &v) {
  double norm = 0;
  for (const auto &value : v) {
    norm += value * value;
  }
  return sqrt(norm);
}

double NormInf(const std::vector<double> &v) {
  double max = std::numeric_limits<double>::lowest();
  for (const auto &value : v) {
    max = std::fabs(value) > max ? std::fabs(value) : max;
  }
  return max;
}

struct Stencil {
  Stencil(double h)
      : C(4.0 / (h * h) + 4 * M_PI * M_PI), N(-1.0 / (h * h)),
        S(-1.0 / (h * h)), W(-1.0 / (h * h)), E(-1.0 / (h * h)) {}
  const double C, N, S, W, E;
};

struct StencilPoint {
  int i;
  int j;
  double weight;
};

enum Cell { UNKNOWN = 0, DIR = 1, NEU = 2, ROB = 0 };

void PoissonJacobiStencil(size_t resolution,
                          size_t iterations) {

  size_t N = resolution + 1;
  double h = 1.0 / resolution;

  const auto stencil = Stencil(h);

  // domain cell types
  std::vector<int> domain(N * N, Cell::UNKNOWN);
  MatrixView<int> domainView(domain, N, N);
  for (size_t n = 0; n != N; ++n) {
    domainView.set(n, 0) = Cell::DIR;
    domainView.set(0, n) = Cell::DIR;
    domainView.set(N - 1, n) = Cell::DIR;
    domainView.set(n, N - 1) = Cell::DIR;
  }

  // referenceSolution
  std::vector<double> referenceSolution(N * N, 0);
  MatrixView<double> referenceSolutionView(referenceSolution, N, N);
  for (size_t j = 0; j != N; ++j) {
    for (size_t i = 0; i != N; ++i) {
      referenceSolutionView.set(i, j) = ParticularSolution(i * h, j * h);
    }
  }

  // right hand side
  std::vector<double> rightHandSide(N * N, 0);
  MatrixView<double> rightHandSideView(rightHandSide, N, N);
  for (size_t j = 0; j != N; ++j) {
    for (size_t i = 0; i != N; ++i) {
      rightHandSideView.set(i, j) =
          ParticularSolution(i * h, j * h) * 4 * M_PI * M_PI;
    }
  }

  auto SolverJacobi = [](std::vector<double> &sol, std::vector<double> &sol2,
                         std::vector<double> &rhs, const Stencil &stencil,
                         size_t N) {
    MatrixView<double> solView(sol, N, N);
    MatrixView<double> sol2View(sol2, N, N);
    MatrixView<double> rhsView(rhs, N, N);

    for (size_t j = 1; j < N - 1; ++j) {
      for (size_t i = 1; i < N - 1; ++i) {
        sol2View.set(i, j) =
            1.0 / stencil.C *
            (rhsView.set(i, j) - (solView.get(i + 1, j) * stencil.E +
                                  solView.get(i - 1, j) * stencil.W +
                                  solView.get(i, j + 1) * stencil.S +
                                  solView.get(i, j - 1) * stencil.N));
      }
    }
    sol.swap(sol2);
  };

  auto ComputeResidual = [](std::vector<double> &sol, std::vector<double> &rhs,
                            const Stencil &stencil, size_t N) {
    MatrixView<double> solView(sol, N, N);
    MatrixView<double> rhsView(rhs, N, N);

    std::vector<double> residual(N * N, 0);
    MatrixView<double> residualView(residual, N, N);
    for (size_t j = 1; j != N - 1; ++j) {
      for (size_t i = 1; i != N - 1; ++i) {
        residualView.set(i, j) =
            rhsView.get(i, j) -
            (solView.get(i, j) * stencil.C + solView.get(i + 1, j) * stencil.E +
             solView.get(i - 1, j) * stencil.W +
             solView.get(i, j - 1) * stencil.S +
             solView.get(i, j + 1) * stencil.N);
      }
    }
    return residual;
  };
  auto ComputeError = [](std::vector<double> &sol,
                         std::vector<double> &reference, size_t N) {
    MatrixView<double> solView(sol, N, N);
    MatrixView<double> referenceView(reference, N, N);

    std::vector<double> error(N * N, 0);
    MatrixView<double> errorView(error, N, N);

    for (size_t j = 1; j != N - 1; ++j) {
      for (size_t i = 1; i != N - 1; ++i) {
        errorView.set(i, j) = referenceView.get(i, j) - solView.get(i, j);
      }
    }
    return error;
  };

  // solution approximation starting with boundary initialized to diriclet
  // conditions, else 0
  std::vector<double> solution(N * N, 0);
  MatrixView<double> solutionView(solution, N, N);
  for (size_t j = 0; j != N; ++j) {
    for (size_t i = 0; i != N; ++i) {
      if (domainView.get(i, j) == Cell::DIR)
        solutionView.set(i, j) = ParticularSolution(i * h, j * h);
      else if (domainView.get(i, j) == Cell::UNKNOWN)
        solutionView.set(i, j) = 0;
      else
        throw("BC not implemented");
    }
  }
  std::vector<double> solution2 = solution;
  std::cout << "solve LSE using stencil jacobi" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t iter = 0; iter <= iterations; ++iter) {
    SolverJacobi(solution, solution2, rightHandSide, stencil, N);
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto seconds =
      std::chrono::duration_cast<std::chrono::duration<double>>(stop - start)
          .count();
  std::cout << std::scientific << "runtime=" << seconds << std::endl;

  {
    auto residual = ComputeResidual(solution, rightHandSide, stencil, N);
    auto residualNorm = NormL2(residual);
    std::cout << std::scientific << "|residual|=" << residualNorm << std::endl;
    auto residualMax = NormInf(residual);
    std::cout << std::scientific << "|residualMax|=" << residualMax
              << std::endl;
    auto error = ComputeError(solution, referenceSolution, N);
    auto errorNorm = NormL2(error);
    std::cout << std::scientific << "|error|=" << errorNorm << std::endl;
    auto errorMax = NormInf(error);
    std::cout << std::scientific << "|errorMax|=" << errorMax << std::endl;
  }
}
