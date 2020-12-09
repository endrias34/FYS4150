
To compile write :
  
```
g++ -fopenmp IsingSocialMain.cpp
```
  
a.out is created  

To run it your first argument will be how big of a system (#n spins)  

./a.out 10 

is 10 spins.

Second argument is if you want the system ordered or not :  

./a.out 10 1   -> is a systems with 10 spins, started with all 1's  
./a.out 10 0   -> is a systems with 10 spins, started -1 or 1 as opinon drawn from unifrom distribution.  

Third argument and fourth argument is to decide if you want to save the system to file at each step inputed    
3 = steps, 4 = magnetization or spins (0 for magnetizations and 1 for spins)
  
./a.out 10 0 10 1  -> 10 spins, unordered, where each spin is saved to file every 10'th MC cycle  
./a.out 10 0 0 0  -> 10 spins, unordered, and now we run 1000 systems and save the magnetizations from all the systems to one file  

Fifth argument is concentration level  

cB * N unordered  :  
./a.out 10 0 10 0 0.9  -> 10 = spins, 0 = unordered, 10 = save every 10'th MC cycle, 0 = Magnetization, 0.9 = 90% spins started with -1  

cB * N with the concentration cB starting in a cluster for index 0 to index cB * N :  
./a.out 10 1 10 0 0.9  -> 10 = spins, 1 = ordered, 10 = save every 10'th MC cycle, 0 = Magnetization, 0.9 = 90% spins started with -1  



sixth argument is chance of not following the rules  
  
./a.out 10 0 10 0 0.9 0.000002  -> 
10 = spins, 0 = unordered, 10 = save every 10'th MC cycle,   
0 = Magnetization, 0.9 = 90% spins started with 1, 0.000002 = 0.00002% chance of not following the rule
  
or with all spins starting with -1  :  
./a.out 10 0 10 0 1 0.000002  -> 
10 = spins, 0 = unordered, 10 = save every 10'th MC cycle,   
0 = Magnetization, 1 = all spins -1, 0.000002 = 0.00002% chance of not following the rule








