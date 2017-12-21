DIR=./
para_sortgenome: para_sortgenome.o
	mpicxx -o para_sortgenome para_sortgenome.o
para_sortgenome.o:para_sortgenome.cc para_sortgenome.h
	mpicxx -c para_sortgenome.cc

.PHONY: clean

clean: 
	rm $(DIR).o para_sortgenome
