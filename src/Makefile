#CC=gcc
CXX=g++
RM=rm -f
#CXXFLAGS=-Wall -O0 -g -fopenmp -pthread -std=c++11 -Wl,--no-as-needed
CXXFLAGS=-Wall -O3 -fopenmp -pthread -lpthread -std=c++11 -Wl,--no-as-needed
LDFLAGS=-fopenmp -pthread -std=c++11 -Wl,--no-as-needed
LDLIBS=-lpthread -lgsl -lblas

SRCS=jacobi.cpp

.PHONY: all clean dist-clean

all: jacobi

jacobi.o: jacobi.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -c $<

jacobi: jacobi.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

clean:
	$(RM) jacobi.o

dist-clean: clean
	$(RM) jacobi 
