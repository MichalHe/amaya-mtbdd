CC=gcc
CXX=g++
CFLAGS=$(shell pkg-config --libs sylvan) -g

amaya-mtbdd.so: mtbdd.cpp amaya_mtbdd.hpp
	$(CXX) -shared -fPIC $(CFLAGS) -o $@ mtbdd.cpp

mtbdd-test: mtbdd.cpp
	$(CXX) $(CFLAGS) -o $@ $?

bdd-test: bdd.c
	$(CC) $(CFLAGS) $? -o $@

