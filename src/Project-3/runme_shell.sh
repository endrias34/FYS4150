#!/bin/bash
for i in case6 # Senario
   do
   for j in verlet_GR  # method
     do
       for k in 100   # Maximum time in yr
          do
            for l in 7 # Time step size 10^-l
               do
                 for m in 1000 # print position every m*10^-l yr
                 do
                    for bb in 2.0 # Inverse square law
                       do
                        ./runme $i $j $k $l $m $bb > Results/Mercury_angle
                    done
                done
            done 
        done
    done
done

# # Full solar system
# for i in case5 # Senario
#  do
#    for j in verlet  # method
#      do
#        for k in 250.0     # Maximum time in yr
#           do
#             for l in 5       # Time step size 10^-l
#                do
#                  for m in 100  # print position every m*10^-l yr
#                    do
#                      for bb in 2.0 # Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb 
#                     done
#                 done
#             done 
#         done
#     done
# done

# #Three body problem
# for i in case4 case41 case411 # #case3 case31 case311  # Senario
#  do
#    for j in verlet  # method
#      do
#        for k in 1000     # Maximum time in yr
#           do
#             for l in 5       # Time step size 10^-l
#                do
#                  for m in 5000  # print position every m*10^-l yr
#                    do
#                      for bb in 2.0 # Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb 
#                     done
#                 done
#             done 
#         done
#     done
# done

# #Escape velocity
# for i in case111d case111c ##  #case111 case111b #  # Senario
#  do
#    for j in verlet  # method
#      do
#        for k in 100.0        # Maximum time in yr
#           do
#             for l in 5       # Time step size 10^-l
#                do
#                  for m in 100  # print position every m*10^-l yr
#                    do
#                      for bb in 3.0 #2.0 # Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb 
#                     done
#                 done
#             done 
#         done
#     done
# done

# # Conservation of energy and angular momentum for elliptical orbit
# for i in case11    # Senario
#  do
#    for j in verletEnergy verletAngmom  # method
#      do
#        for k in 0.5        # Maximum time in yr
#           do
#             for l in 5       # Time step size 10^-l
#                do
#                  for m in 10  # print position every m*10^-l yr
#                    do
#                      for bb in 3.0 #2.0  Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb > Results/Ellipse_$bb"_"$j.txt
#                     done
#                 done
#             done 
#         done
#     done
# done
# rm Results/earth_*Energy*
# rm Results/earth_*Angmom*

# # Modified gravity
# for i in case11 #case1 #    # Senario
#  do
#    for j in euler verlet    # method
#      do
#        for k in 10.0        # Maximum time in yr
#           do
#             for l in 5       # Time step size 10^-l
#                do
#                  for m in 100  # print position every m data point
#                    do
#                      for bb in 2.0  2.01 2.5 2.95 3.0 # Gravitational law beta
#                         do
#                         ./runme $i $j $k $l $m $bb
#                     done
#                 done
#             done 
#         done
#     done
# done

# # CPU time for forward Euler and velocity Verlet Algorithms
# # ./runme_shell.sh > Results/eulerTime.txt
# # ./runme_shell.sh > Results/verletTime.txt
# for i in case1    # Senario
#  do
#    for j in verletTime #eulerTime # # method
#      do
#        for k in 1.0        # Maximum time in yr
#           do
#             for l in 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7 2 3 4 5 6 7    # Time step size 10^-l
#                do
#                  for m in 100000  # print position every m*10^-l yr
#                    do
#                      for bb in 2.0  # Inverse square law
#                         do

#                         ./runme $i $j $k $l $m $bb 

#                     done
#                 done
#             done 
#         done
#     done
# done
# rm Results/earth_*Time*

# # Conservation of energy and angular momentum
# for i in case1    # Senario
#  do
#    for j in eulerEnergy verletEnergy eulerAngmom verletAngmom  # method
#      do
#        for k in 50.0        # Maximum time in yr
#           do
#             for l in 5       # Time step size 10^-l
#                do
#                  for m in 100  # print position every m*10^-l yr
#                    do
#                      for bb in 2.0  # Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb > Results/$j.txt
#                     done
#                 done
#             done 
#         done
#     done
# done
# rm Results/earth_*Energy*
# rm Results/earth_*Angmom*

# # Circular orbit accuracy
# for i in case1    # Senario
#  do
#    for j in euler verlet    # method
#      do
#        for k in 1.0        # Maximum time in yr
#           do
#             for l in 2 3 4 5 6 7      # Time step size 10^-l
#                do
#                  for m in 10  # print position every m*10^-l yr
#                    do
#                      for bb in 2.0  # Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb
#                     done
#                 done
#             done 
#         done
#     done
# done

# # Circular orbit test
# for i in case1    # Senario
#  do
#    for j in euler verlet    # method
#      do
#        for k in 50.0        # Maximum time in yr
#           do
#             for l in 3       # Time step size 10^-l
#                do
#                  for m in 10  # print position every m*10^-l yr
#                    do
#                      for bb in 2.0  # Inverse square law
#                         do
#                         ./runme $i $j $k $l $m $bb
#                     done
#                 done
#             done 
#         done
#     done
# done