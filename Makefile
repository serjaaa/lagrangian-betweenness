#To copy this file on Niko
#scp -P 8022 ~/cpp/code_Betweenness/Makefile abaudena@localhost:~/cpp/code_Betweenness/

#CC_MONO=g++ -Wall
#CC_CHECK=g++ -Wall -C -E
CC_SIMPLE=g++

#CC=$(CC_MONO)
#CC=$(CC_CHECK)
CC=${CC_SIMPLE}

LIBS=-lnetcdf_c++
RM='rm -rf'

all: betw_computation.out

parameters.o: parameters.cpp parameters.h
	${CC} -c parameters.cpp

read_velocity.o: read_velocity.cpp read_velocity.h
	${CC} -c read_velocity.cpp ${LIBS}

betw_computation.o: betw_computation.cpp parameters.h read_velocity.h
	${CC} -c betw_computation.cpp -fopenmp

betw_computation.out: betw_computation.o parameters.o read_velocity.o parameters.h read_velocity.h
	${CC} betw_computation.o parameters.o read_velocity.o -o betw_computation.out ${LIBS} -fopenmp

clean:
	${RM} *.o 
