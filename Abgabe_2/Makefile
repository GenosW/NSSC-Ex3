# requirements on ubuntu/linux mint

CXX=g++
# Single thread -> no omp
#CXXFLAGS=-std=c++11 -Wall -pedantic -O3 -ffast-math -march=native
# Linux
CXXFLAGS= -v -std=c++14 -O3 -Wall -pedantic -march=native -ffast-math -fopenmp
# Mac
#CXXFLAGS= -v -std=c++14 -O3 -Wall -pedantic -march=native -ffast-math -Xpreprocessor -fopenmp -lomp  

.DEFAULT_GOAL := all

all: cg
 
cg: Task_2.cpp Makefile
	$(CXX) Task_2.cpp -o cg $(CXXFLAGS)

.PHONY: clean
clean:
	rm cg
