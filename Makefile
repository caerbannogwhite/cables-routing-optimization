P=wfsolver
SHELL=/bin/sh

CC=g++
CFLAGS=-g -Wall -O2
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128/cplex
LDLIBS=-L$(CPLEX_HOME)/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -ldl
INC=-I$(CPLEX_HOME)/include/ilcplex

$(P) : aux.o heur.o main.o wf.o
	$(CC) $(CFLAGS) bin/aux.o bin/heur.o bin/main.o bin/wf.o -o $@ $(LDLIBS) $(INC)
	mv $@ bin

aux.o : src/$*.cpp src/aux.h src/wf.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

heur.o : src/$*.cpp src/heur.h src/aux.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

main.o : src/$*.cpp src/wf.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

wf.o : src/$*.cpp src/wf.h
	$(CC) $(CFLAGS) -c src/$*.cpp $(LDLIBS) $(INC)
	mv $@ bin

clean:
	rm -f bin/$(P) bin/*.o

.PHONY: clean