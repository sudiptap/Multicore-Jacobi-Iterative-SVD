#CC=gcc
CXX=g++
RM=rm -f
#CXXFLAGS=-Wall -O0 -g -fopenmp -pthread -std=c++11 -Wl,--no-as-needed
CXXFLAGS=-Wall -O3 -fopenmp -pthread -lpthread -std=c++11 -Wl,--no-as-needed
LDFLAGS=-fopenmp -pthread -std=c++11 -Wl,--no-as-needed
LDLIBS=-lpthread -lgsl -lblas

SRCS=jacobi_svd.cpp jacobi_async.hpp jacobi_gsl.hpp jacobi_async_delayed.hpp jacobi_gsl_random_pair.hpp jacobi_gsl_random_permutation.hpp jacobi_gsl_optimal_pair.hpp jacobi_gsl_best_pair_first.hpp

.PHONY: all clean dist-clean

all: jacobi jmos

jacobi_svd.o: jacobi_svd.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -c $<

jacobi: jacobi_svd.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

jmos.o: jacobi_matrix_one_sided.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

jmos: jmos.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

clean:
	$(RM) jacobi_svd.o jmos.o

dist-clean: clean
	$(RM) jacobi jmos
