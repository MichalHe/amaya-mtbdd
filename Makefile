CXX=g++
CXXLIBS=$(shell pkg-config --libs sylvan)
CXXFLAGS=-O2 -shared -fPIC

.PHONY := clean

shared-lib: build build/amaya-mtbdd.so

build/amaya-mtbdd.so: build/wrapper.o build/operations.o build/custom_leaf.o build/base.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CXXLIBS)

build/wrapper.o: src/wrapper.cpp include/wrapper.hpp include/operations.hpp include/base.hpp include/custom_leaf.hpp 
	$(CXX) -c $(CXXFLAGS) src/wrapper.cpp -o build/wrapper.o

build/operations.o: src/operations.cpp include/operations.hpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/operations.cpp -o build/operations.o

build/custom_leaf.o: src/custom_leaf.cpp include/custom_leaf.hpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/custom_leaf.cpp -o build/custom_leaf.o

build/base.o: src/base.cpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/base.cpp -o build/base.o

build:
	-mkdir build

clean:
	rm base.o custom_leaf.o operations.o wrapper.o amaya-mtbdd.so
