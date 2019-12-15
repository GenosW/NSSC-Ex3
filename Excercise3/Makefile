CXX=g++
# For Linux
CXXFLAGS=-std=c++11 -Wall -pedantic -O3 -march=native -fopenmp
  
.DEFAULT_GOAL := all

all: mcint 
 
mcint: mcint.cpp Makefile
	$(CXX) mcint.cpp -o mcint -v  $(CXXFLAGS)  

.PHONY: clean
clean:
	rm mcint