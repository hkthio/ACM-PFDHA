# Makefile for the AECOM PFDHA code base

FFLAGS= -O3 -g -Mextend
CFLAGS = -O3 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/opt/local/include

F77= gfortran
CC = /usr/local/bin/gcc

.c.o: 
	${CC}  $(CFLAGS) -c $*.c

.f.o: 
	${F77}  $(FFLAGS) -c $*.f

PHA-IAEA: PHA-IAEA.o PHA-Disp.o PHA-Math.o
	${CC} -o PHA-IAEA PHA-IAEA.o PHA-Disp.o PHA-Math.o

PFDHA-IAEA: PFDHA-IAEA.o PFDHA-Disp.o PHA-Math.o
	${CC} -o PFDHA-IAEA PFDHA-IAEA.o PFDHA-Disp.o PHA-Math.o

PHA-Test: PHA-Test.o PHA-Disp.o PHA-Math.o
	${CC} -o PHA-Test PHA-Test.o PHA-Disp.o PHA-Math.o

PHA-HazCurve: PHA-HazCurve.o PHA-Disp.o PHA-Math.o
	${CC} -o PHA-HazCurve PHA-HazCurve.o PHA-Disp.o PHA-Math.o

PHA-TestD: PHA-TestD.o PHA-Disp.o PHA-Math.o
	${CC} -o PHA-TestD PHA-TestD.o PHA-Disp.o PHA-Math.o

