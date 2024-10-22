parallel: parallel.cu cuda_cell.h parallel.h
	-mkdir exe
	nvcc parallel.cu -o parallel.o -dc
	nvcc parallel.o -o exe/parallel -lm
	-rm parallel.o

run-parallel: exe/parallel
	./exe/parallel

sequential: sequential.c sequential/functions.c cell.h sequential/functions.h
	-mkdir exe
	gcc sequential.c sequential/functions.c sequential/simulation_utils.c sequential/physics.c -o exe/sequential -lm

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