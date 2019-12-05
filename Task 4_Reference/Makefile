# requirements on ubuntu
# sudo apt-get build-essentials
CXX=g++

# Update the path to your Eigen library
EIGEN=/iue/tcad3x/users/XXX

# Single thread -> no omp
#CXXFLAGS=-std=c++11 -Wall -pedantic -O3 -ffast-math -march=native
# Linux
#CXXFLAGS=-std=c++14 -O3 -Wall -pedantic -march=native -ffast-math -fopenmp
# Mac
CXXFLAGS= -v -std=c++14 -O3 -Wall -pedantic -march=native -ffast-math #-Xpreprocessor -fopenmp -lomp  

.DEFAULT_GOAL := all

all: stenciljacobi
 
stenciljacobi: Makefile solver.cpp poisson_jacobi.hpp arguments.hpp
	$(CXX) solver.cpp -o stenciljacobi -lpthread -v $(CXXFLAGS) -I$(EIGEN)

.PHONY: clean
clean:
	rm stenciljacobi
