CXX=g++
CXXLIBS=$(shell pkg-config --libs sylvan)

SHARED_LIB_FLAGS=-shared -fPIC
CXXFLAGS=-O2 -g --std=c++20 $(SHARED_LIB_FLAGS)

.PHONY := clean

lazy-tests: build/lazy-tests.o build/lazy.o build/base.o build/custom_leaf.o build/operations.o
	$(CXX) -o $@ $^ $(CXXLIBS)

build/lazy-tests.o: src/lazy-tests.cpp
	$(CXX) -c $(CXXFLAGS) -I ./external -o $@ $^

shared-lib: build build/amaya-mtbdd.so

build/lazy.o: build src/lazy.cpp include/lazy.hpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/lazy.cpp -o build/lazy.o

build/amaya-mtbdd.so: build/wrapper.o build/operations.o build/custom_leaf.o build/hopcroft_leaf.o build/base.o build/lazy.o
	$(CXX) $(CXXFLAGS) $(SHARED_LIB_FLAGS) -o $@ $^ $(CXXLIBS)

build/wrapper.o: src/wrapper.cpp include/wrapper.hpp include/operations.hpp include/base.hpp include/custom_leaf.hpp include/hopcroft_leaf.hpp
	$(CXX) -c $(CXXFLAGS) src/wrapper.cpp -o build/wrapper.o

build/operations.o: src/operations.cpp include/operations.hpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/operations.cpp -o build/operations.o

build/custom_leaf.o: src/custom_leaf.cpp include/custom_leaf.hpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/custom_leaf.cpp -o build/custom_leaf.o

build/hopcroft_leaf.o: src/hopcroft_leaf.cpp include/hopcroft_leaf.hpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/hopcroft_leaf.cpp -o build/hopcroft_leaf.o

build/base.o: src/base.cpp include/base.hpp
	$(CXX) -c $(CXXFLAGS) src/base.cpp -o build/base.o

build:
	-mkdir build

clean:
	rm build/* || true

build/test.o: src/test.cpp include/hopcroft_leaf.hpp
	$(CXX) -c $(CXXFLAGS) src/test.cpp -o build/test.o

test: build/test.o build/hopcroft_leaf.o build/wrapper.o build/operations.o build/custom_leaf.o build/base.o
	$(CXX) -o $@ $^ $(CXXLIBS)
