

To run the program put all files in a folder and write  
  
make  
  
on the commandline in your shell.

Makefile use g++ version for compiling and if you do not have it, open Makefile change line 4  
to something like  
  
CPPflags= c++ -fopenmp -O3 -std=c++11  
  
After compiling write ./runme in the commandline and you will get alternatives where one is you  
setting your own configurations, and the other one is the configurations set between line 71 and  
93 in mainprogram.cpp.  


If you want to configure through commandline write  

./runme 1  

and you will get an error where what to be put in is give  

