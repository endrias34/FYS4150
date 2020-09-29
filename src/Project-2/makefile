# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas

# -------------------------------------------------------------------------------- #
# Here we define the name of the executable
# PROG = runme_beam
# ${PROG} :	initialize.o  main_beam.o  jacobi.o
# 			${CPPflags} initialize.o main_beam.o jacobi.o ${LIB} -o ${PROG}

# initialize.o :		initialize.cpp
# 		        	${CPPflags} -c initialize.cpp

# main_beam.o :		main_beam.cpp 
# 		        	${CPPflags} -c main_beam.cpp

# jacobi.o :		    jacobi.cpp 
# 		        	${CPPflags} -c jacobi.cpp


# -------------------------------------------------------------------------------- #

# Here we define the name of the executable
# PROG = runme_one_electron
# ${PROG} :	initialize.o  main_one_electron.o  jacobi.o
# 			${CPPflags} initialize.o main_one_electron.o jacobi.o ${LIB} -o ${PROG}

# initialize.o :		initialize.cpp
# 		        	${CPPflags} -c initialize.cpp

# main_one_electron.o :		main_one_electron.cpp 
# 		        	${CPPflags} -c main_one_electron.cpp

# jacobi.o :		    jacobi.cpp 
# 		        	${CPPflags} -c jacobi.cpp

# -------------------------------------------------------------------------------- #

# Here we define the name of the executable
PROG = runme_two_electrons
${PROG} :	initialize.o  main_two_electrons.o  jacobi.o
			${CPPflags} initialize.o main_two_electrons.o jacobi.o ${LIB} -o ${PROG}

initialize.o :		initialize.cpp
		        	${CPPflags} -c initialize.cpp

main_two_electrons.o :		main_two_electrons.cpp 
		        	${CPPflags} -c main_two_electrons.cpp

jacobi.o :		    jacobi.cpp 
		        	${CPPflags} -c jacobi.cpp

# -------------------------------------------------------------------------------- #

#PROG = runme_test
#${PROG} :	unit_tests.o  test_main.o  jacobi.o
#			${CPPflags} unit_tests.o test_main.o jacobi.o ${LIB} -o ${PROG}

#unit_tests.o :		unit_tests.cpp
#		        	${CPPflags} -c unit_tests.cpp

#test_main.o :		test_main.cpp 
#		        	${CPPflags} -c test_main.cpp

#jacobi.o :		    jacobi.cpp 
#		        	${CPPflags} -c jacobi.cpp

