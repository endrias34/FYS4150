# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3 -std=c++11
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas

# -------------------------------------------------------------------------------- #
# Here we define the name of the executable
PROG = runme
${PROG} :	main.o  planet.o  solver.o
			${CPPflags} main.o planet.o	solver.o ${LIB} -o ${PROG}

main.o	:	main.cpp 
			${CPPflags}	-c	main.cpp

planet.o :		 planet.cpp
					${CPPflags}	-c	planet.cpp

solver.o :  solver.cpp
					${CPPflags}	-c	solver.cpp

clean:
	rm *.o