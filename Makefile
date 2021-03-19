CC=gcc
CXX=g++
CFLAGS=$(shell pkg-config --libs sylvan) -g

mtbdd-test: mtbdd.cpp
	$(CXX) $(CFLAGS) -o $@ $?

bdd-test: bdd.c
	$(CC) $(CFLAGS) $? -o $@

