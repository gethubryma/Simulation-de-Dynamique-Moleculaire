CXX = g++
CXXFLAGS = -Wall -O3 -std=c++11 -fopenmp

all: simulation

simulation: main.cpp Simulation.cpp Simulation.h
	$(CXX) $(CXXFLAGS) main.cpp Simulation.cpp -o simulation

clean:
	rm -f simulation
