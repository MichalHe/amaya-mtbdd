CC=gcc
CXX=g++
CFLAGS=$(shell pkg-config --libs sylvan) -O2 -shared -fPIC

amaya-mtbdd.so: wrapper.o operations.o custom_leaf.o base.o
	$(CXX) -o $@ $^ $(CFLAGS)

wrapper.o: wrapper.cpp wrapper.hpp operations.hpp base.hpp custom_leaf.hpp 
	$(CXX) -c wrapper.cpp $(CFLAGS) -o wrapper.o

operations.o: operations.cpp operations.hpp base.hpp
	$(CXX) -c operations.cpp $(CFLAGS) -o operations.o

custom_leaf.o: custom_leaf.cpp custom_leaf.hpp base.hpp
	$(CXX) -c custom_leaf.cpp $(CFLAGS) -o custom_leaf.o

base.o: base.cpp base.hpp
	$(CXX) -c base.cpp $(CFLAGS) -o base.o

clean:
	rm base.o custom_leaf.o operations.o wrapper.o amaya-mtbdd.so
