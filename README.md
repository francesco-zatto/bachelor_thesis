# HIS simulation in CUDA

In the repo there are different folders with C and CUDA source codes to run simulations, tests of the HIS, Human Immune System. 

Thanks to the main Makefile it is possible to compile and test each available source code. Of course a change in a line can cause a change in the Makefile if it requires a dependency of other files.

## Langevin equation: langevin folder
First thing first to test for an HIS simulation is the movement of a cell or a molecule. The movement can be modelled following Langevin equations to create a brownian movement. It creates a random pattern, because a particle in a fluid collides non-stop with every particle nearby.
There are two ways to draw the brownian movement:

- using continous functions, where in each timestep a particle can take any coordinate in the space. In this case, to try Langevin continous equation, first it is mandatory to compile the C source code and then it is possible to execute the resulting executable.
 ```
 make langevin          #to compile
 make run-langevin      #to run 
 ```

 - using discrete functions, where in each timestep a particle can take a coordinate with integer values, as a place inside a space grid, like the final HIS CUDA simulation. As mentioned before, it has to be compiled and then executed.
  ```
 make discrete-langevin          #to compile
 make discrete-run-langevin      #to run 

 #to run using a certain value for the number of cells to move
 ./langevin/exe/discrete_langevin <cells>
 ```

 In both cases, the final result is a `csv` file containing the particle position taken in each timestep. Then you can use it as you like or you can create a dataframe and then plot the result using the jupyter notebook file `plotter_cell.ipynb`. There is still a lot of manual work to do, it could be improved taking automatically the csv result and plot it in an image.

 ## Tests: tests folder

 At the moment the only test written and executed is for the random generation of the cells inside the discrete 2D grid.
```
 make test-generation          #to compile
 make run-test-generation      #to run with default values

 #to run using certain values for the quantities of B cells, T cells and antigens
 ./tests/exe/generation <B_cells> <T_cells> <Ag>
 ```
 
 Even in this case the result is the plot of the generated cells inside the grid as a `csv` file. And you have to use the `ipynb` plotter file to represent it in a plot.

 ## Simulation: repo folder
 In the main folder of the repository, you can run the sequential or the parallel code.
 In both cases the user can choose the default values used in simulation, the ones present in `simulation_utils.h`, or they can choose any value they want for grid's size or the number of cells. The next commands can be ran only when the required executable is already present. To check how to compile the source codes check the next section for sequential and parallel programs. 
 ```
 # To run it with default values
 make run-sequential
 make run-parallel

 #To run it with a chosen grid size
 ./exe/sequential <n>
 ./exe/parallel <n>

 #To run it with a chosen number of B cells, T cells and antigens
 ./exe/sequential <b> <t> <ag>
 ./exe/parallel <b> <t> <ag>

 #To run it with a chosen grid size and chosen number of cells
 ./exe/sequential <b> <t> <ag> <n>
 ./exe/parallel <b> <t> <ag> <n>
 ```

### Sequential program
In the sequential program the whole simulation will be ran without any piece of parallel code, so no CUDA at all. To compile it just run the following comand:
```
make sequential
```

### Parallel program
For the parallel implementation the nvcc compiler is required and only in that case it is possible to compile it using the comand:
```
make parallel
```

### Timing tests
Inside the main folder of the repository, there are two bash scripts, `sequential_tests.sh` and `parallel_tests.sh`, that compile and run different tests using different input parameters of the two versions. The results of the tests, only the execution time of the simulation, is saved in two csv files inside the results folder, `seq_tests.csv` and `parallel_tests.csv`, and then they're used to compare and plot timings with `plot_results.ipynb`.