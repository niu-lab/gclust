#FLAGS = -I./ -O3 -pg
FLAGS = -I ./ -O3
SRC = gclust.cpp paraSA.cpp fasta.cpp

all: gclust 

gclust: gclust.o paraSA.o fasta.o
	g++   $(FLAGS) $^ -o $@ -lpthread

.cpp.o:
	g++   $(FLAGS) -Wall -c $<

.c.o:
	gcc $(FLAGS) -Wall -c $<

clean: 
	rm -f *.o gclust

