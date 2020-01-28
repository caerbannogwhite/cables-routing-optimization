P=wfsolver
SHELL=/bin/sh

CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128/cplex
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2

LDLIBS=-L$(CPLEX_HOME)/lib/x86-64_linux/static_pic -lcplex -L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl
INC=-I$(CPLEX_HOME)/include/ilcplex -I$(BOOST_HOME)

$(P) : aux.o heur.o main.o wf.o
	$(CC) $(CFLAGS) bin/aux.o bin/heur.o bin/main.o bin/wf.o -o $@ $(LDLIBS) $(INC)
	mv $@ bin

aux.o : src/aux.cpp src/aux.h src/wf.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

heur.o : src/heur.cpp src/heur.h src/aux.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

main.o : src/main.cpp src/wf.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

wf.o : src/wf.cpp src/wf.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

clean:
	rm -f bin/$(P) bin/*.o

.PHONY: clean