The C++ file is executed after compilation by "./executionfile 10", where 10 in this case is your choice of gridpoints. 

The python file is ready for plotting the output files created by running you execution file for 10, 100, 1000 as input. If you want to experiment with other stepsizes, be sure to edit line 7

```
problem_sizes = np.logspace(1, 3, 3, endpoint = True).astype(int)
```

to include your the values you tried when you want to plot it.
If you want to see the full effect of where the precision starts to degrad, because of computational limits, comment out line 46 in Poission1D.cpp

```
sol_lu(n);   // LU-decomposition based algorithm
```

make an new executable file, and run it with from 10^n, n >= 5, and higher.
After you have tried your desired values, plot their realtive error by uncommenting line 92 in plot.py

```
# problem_sizes = np.logspace(1, n, n, endpoint = True).astype(int)
```

and inserting the maximum n you have tried. If you have not gone step by step in your n's, do it like this

```
problem_sizes = [10, 100, 1000, valuetried1, valuetried2, ...]
```
