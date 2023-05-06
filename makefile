CXX ?= g++
CXXFLAGS = -std=c++17 -Wall
OPTFLAGS ?= -O3    -flto  -fomit-frame-pointer  -march=native 
OPENMPFLAG = -fopenmp

all: kakuro_solver

kakuro_solver: kakuro_solver.o
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(OPENMPFLAG) kakuro_solver.o -o kakuro_solver

kakuro_solver.o: kakuro_solver.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(OPENMPFLAG) -c kakuro_solver.cpp

clean:
	rm -f *.o kakuro_solver