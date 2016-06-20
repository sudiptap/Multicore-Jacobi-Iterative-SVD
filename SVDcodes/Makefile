#CC=gcc
CXX=g++
RM=rm -f
#CXXFLAGS=-Wall -O0 -g -fopenmp -pthread -std=c++11 -Wl,--no-as-needed
#CXXFLAGS=-Wall -O3 -fopenmp -pthread -lpthread -std=c++11 -Wl,--no-as-needed
#CXXFLAGS=-Wall -O3 -std=c++11 
CXXFLAGS=-Wall -O0 -g -std=c++11 
LDFLAGS=-fopenmp -pthread -std=c++11 -Wl,--no-as-needed
LDLIBS=-lpthread -lgsl -lblas

SRCS=svd.h \
		 cyclic_jacobi_one.cpp \
		 cyclic_jacobi_two.cpp \
		 jrs_one.cpp \
		 jrs_two.cpp \
		 jps_one.cpp \
		 jps_two.cpp \
		 jprs_one.cpp \
		 jprs_two.cpp \
		 par_jrs_one.cpp \
		 par_jrs_two.cpp \
		 par_jps_one.cpp \
		 par_jps_two.cpp \
		 par_jprs_one.cpp \
		 par_jprs_two.cpp \


#jacobi_gsl_best_pair_first.hpp \

.PHONY: all clean

all: svd 

svd: svd.cpp ${SRCS}
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS) 

clean:
	$(RM) svd
