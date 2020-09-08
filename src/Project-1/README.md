The C++ file is executed after compilation by "./executionfile 10", where 10 in this case is your choice of gridpoints. 

The python file is set ready for the output file created by running you execution file for 10, 100, 1000 as input. If you want to experiment with other stepsizes, be sure to edit 
```
problem_sizes = np.logspace(1, 3, 3, endpoint = True).astype(int)
```
to include your desired values when you want to plot it.
