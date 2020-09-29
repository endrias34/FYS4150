#!/bin/bash

# ------- Run Beam ----------
#for i in $(seq 10 10 250) # Loop over N
#      ./runme_beam $i
#   do
#done

# ------- One Electron ----------
#for i in $(seq 10 10 250) # Loop over N
#   do
#    for j in $(seq 1 1 6) # Loop over rho_max
#     do
#       ./runme_one_electron $i $j 
#   done
#done

# ------- Two Electron ----------
#for i in $(seq 10 10 250) # Loop over N
#   do
#    for j in $(seq 1 1 6) # Loop over rho_max
#     do
#        for k in 0.01 0.05 0.25 0.5 0.01827 1 5 # Loop over omega_r
#          do
#          ./runme_two_electrons $i $j $k 
#        done
#    done
#done
