# COMPHY: path of COMPHY directory
COMPHY	=	$(HOME)/comphy
FC	=	gfortran
#LFLAGS 	= 	-O -L$(COMPHY)/lib -lnumer
LFLAGS 	= 	-g -Llib -lnumer
FFLAGS 	= 	-c -g

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) 
CXXFLAGS  += $(ROOTCFLAGS)
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3
GXXd	   = /usr/bin/g++ -Wall -g

ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
P5640FLAGS  = -L${P5640LIB}/lib -lP5640  -I${P5640LIB}
GSLFLAGS    = -I${EBROOTGSL}/include/gsl  -I/usr/include/gsl -lgsl -lgslcblas




all: LaplaceLine part1 LaplaceLine1 LaplaceLine2 LaplaceLine3


LaplaceLine: LaplaceLine.C
	g++ -g -Wall -o LaplaceLine LaplaceLine.C $(ROOTFLAGS) $(GSLFLAGS)

part1: part1.C
	g++ -g -Wall -o part1 part1.C $(ROOTFLAGS) $(GSLFLAGS)

part2: part2.C
	g++ -g -Wall -o part2 part2.C $(ROOTFLAGS) $(GSLFLAGS)

LaplaceLine1: LaplaceLine1.C
	g++ -g -Wall -o LaplaceLine1 LaplaceLine1.C $(ROOTFLAGS) $(GSLFLAGS)

LaplaceLine2: LaplaceLine2.C
	g++ -g -Wall -o LaplaceLine2 LaplaceLine2.C $(ROOTFLAGS) $(GSLFLAGS)

LaplaceLine3: LaplaceLine3.C
	g++ -g -Wall -o LaplaceLine3 LaplaceLine3.C $(ROOTFLAGS) $(GSLFLAGS)

clean:
	rm -f LaplaceLine part1 part2 LaplaceLine1 LaplaceLine2 LaplaceLine3 *.o *.so *.pcm *.d *~
