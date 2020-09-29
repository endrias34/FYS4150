#!/bin/bash
#for i in $(seq 10 10 250)
for i in 250
   do

   #./runme_beam $i 

   #for j in $(seq 1 1 6)
   for j in 5 10 20 30 40 50
   #for j in 20
     do

       #./runme_one_electron $i $j 

    #   for k in 0.01 0.5 1 5
        for k in 0.01 0.05 0.25 0.5 0.01827 1 5
          do

          ./runme_two_electrons $i $j $k 

        done
    done
done