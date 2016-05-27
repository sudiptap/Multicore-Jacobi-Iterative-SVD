#CC=gcc
CXX=g++
RM=rm -f
CXXFLAGS=-Wall -O0 -g -fopenmp -pthread -std=c++11 -Wl,--no-as-needed
#CXXFLAGS=-Wall -O3 -fopenmp -pthread -lpthread -std=c++11 -Wl,--no-as-needed
LDFLAGS=-fopenmp -pthread -std=c++11 -Wl,--no-as-needed
LDLIBS=-lpthread -lgsl -lblas

SRCS=jacobi.cpp jacobi_async.cpp jacobi_gsl.cpp

.PHONY: all clean dist-clean

all: jacobi jmos

jacobi.o: jacobi.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -c $<

jacobi: jacobi.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

jmos.o: jacobi_matrix_one_sided.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

jmos: jmos.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

clean:
	$(RM) jacobi.o jmos.o

dist-clean: clean
	$(RM) jacobi jmos
