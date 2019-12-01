# requirements on ubuntu/linux mint
# sudo apt-get install build-essential libopenblas-base libopenblas-dev libeigen3-dev

CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic -O3 -ffast-math -march=native  

.DEFAULT_GOAL := all

all: stenciljacobi
 
stenciljacobi: stenciljacobi.cpp Makefile
	$(CXX) stenciljacobi.cpp -o stenciljacobi -v  $(CXXFLAGS)

.PHONY: clean
clean:
	rm stenciljacobi