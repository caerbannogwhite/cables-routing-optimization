P=wfsolver
SHELL=/bin/sh

CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128/cplex
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2

LDLIBS=-L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl
LDLIBSCPX=-L$(CPLEX_HOME)/lib/x86-64_linux/static_pic -lcplex $(LDLIBS)
INC=-I$(CPLEX_HOME)/include/ilcplex -I$(BOOST_HOME)

$(P) : $(P)_heur $(P)_cpx
	# fake rule

$(P)_heur : common.o main_heur.o wf_heur.o
	$(CC) $(CFLAGS) bin/common.o bin/main_heur.o bin/wf_heur.o -o $@ $(LDLIBS) $(INC)
	mv $@ bin

common.o : src/common.cpp src/common.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

main_heur.o : src/main_heur.cpp src/wf_heur.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

wf_heur.o : src/wf_heur.cpp src/wf_heur.hpp src/common.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

$(P)_cpx : common.o main_cpx.o wf_cpx.o
	$(CC) $(CFLAGS) bin/common.o bin/main_cpx.o bin/wf_cpx.o -o $@ $(LDLIBSCPX) $(INC)
	mv $@ bin

main_cpx.o : src/main_cpx.cpp src/wf_cpx.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBSCPX) $(INC)
	mv $@ bin

wf_cpx.o : src/wf_cpx.cpp src/wf_cpx.hpp src/common.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBSCPX) $(INC)
	mv $@ bin

clean:
	rm -f bin/$(P)_cpx bin/$(P)_heur bin/*.o

.PHONY: clean