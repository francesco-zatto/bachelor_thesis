sequential: sequential.c functions.c cell.h functions.h
	gcc sequential.c -o exe/sequential

run-sequential: exe/sequential
	./exe/sequential

langevin: langevin/langevin.c cell.h
	gcc langevin/langevin.c -o langevin/exe/langevin -lm

run-langevin: langevin/exe/langevin
	./langevin/exe/langevin

discrete-langevin: langevin/discrete_langevin.c cell.h
	gcc langevin/discrete_langevin.c -o langevin/exe/discrete_langevin -lm

run-discrete-langevin: langevin/exe/discrete_langevin
	./langevin/exe/discrete_langevin

test-generation: tests/generation.c 
	gcc tests/generation.c -o tests/exe/generation

run-test-generation: tests/exe/generation
	./tests/exe/generation 

clear:
	-rm langevin/exe/* tests/exe/* exe/* *.csv tests/*.csv langevin/*.csv