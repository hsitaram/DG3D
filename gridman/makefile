FILES=Cell.o Gridman.o Face.o Legendre.o Lagrange.o

grdriver:	${FILES}
	${CXX} ${FILES} grid_driver.cpp -o grdriver.exe

clean:
	rm -rf *.o *.exe *.log *.vtu
