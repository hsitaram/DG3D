GRIDMANDIR=../gridman
CONVDIFFDIR = ../Acoustics
BUILDDIR = .


WFLAGS = -Wall
WFLAGS =
CXXFLAGS = -I${GRIDMANDIR} ${WFLAGS}

GRIDFILES = ${GRIDMANDIR}/Cell.o ${GRIDMANDIR}/Gridman.o ${GRIDMANDIR}/Face.o ${GRIDMANDIR}/Legendre.o ${GRIDMANDIR}/Lagrange.o
CONVDIFFFILES = ${CONVDIFFDIR}/Euler_manager.o ${CONVDIFFDIR}/Euler_timestepper.o ${CONVDIFFDIR}/Euler_outputs.o

BUILDOBJFILES = ${BUILDDIR}/Cell.o ${BUILDDIR}/Gridman.o ${BUILDDIR}/Face.o ${BUILDDIR}/Legendre.o ${BUILDDIR}/Lagrange.o \
 ${BUILDDIR}/Euler_manager.o ${BUILDDIR}/Euler_timestepper.o ${BUILDDIR}/Euler_outputs.o

convdiff:	${GRIDFILES} ${CONVDIFFFILES}
	mv ${GRIDFILES} ${BUILDDIR}
	mv ${CONVDIFFFILES} ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${BUILDOBJFILES} ${CONVDIFFDIR}/Euler_main.cpp -o euler.exe

clean:
	rm -rf *.o *.exe
