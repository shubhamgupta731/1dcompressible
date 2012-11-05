OBJS = Cell.o FreeFunctions.o RiemannSolver.o FluxComputation.o Godunov.o main.o
CC = g++
DEBUG = -g
INCLUDE = -I /usr/local/include/eigen3/ 
CFLAGS = -Wall -c $(DEBUG) $(INCLUDE)
LFLAGS = -Wall $(DEBUG)

Riemann : $(OBJS)
	    $(CC) $(LFLAGS) $(OBJS) -o Riemann

Cell.o : Cell.h Cell.cpp
	    $(CC) $(CFLAGS) Cell.cpp

FreeFunctions.o : Cell.h FreeFunctions.h FreeFunctions.cpp
		$(CC) $(CFLAGS) FreeFunctions.cpp

FluxComputation.o : Cell.h RiemannSolver.h FluxComputation.h FluxComputation.cpp
	    $(CC) $(CFLAGS) FluxComputation.cpp

Godunov.o : FluxComputation.h Godunov.h Godunov.cpp 
	    $(CC) $(CFLAGS) Godunov.cpp

RiemannSolver.o : Cell.h FreeFunctions.h RiemannSolver.h RiemannSolver.cpp
	    $(CC) $(CFLAGS) RiemannSolver.cpp

#HLLC.o : RiemannSolver.h FluxComputation.h HLLC.h HLLC.cpp 
#	    $(CC) $(CFLAGS) HLLC.cpp

main.o : Cell.o Godunov.h HLLC.h FreeFunctions.h main.cpp 
		$(CC) $(CFLAGS) main.cpp

clean:
	    \rm *.o *~ Riemann

