CXX = g++
CXXFLAGS = -std=c++11 -Wall -fopenmp

all: kakuro_solver

kakuro_solver: kakuro_solver.o
	$(CXX) $(CXXFLAGS) kakuro_solver.o -o kakuro_solver

kakuro_solver.o: kakuro_solver.cpp
	$(CXX) $(CXXFLAGS) -c kakuro_solver.cpp

clean:
	rm -f *.o kakuro_solver
