# Tester
Simulation project for modeling metacommunity dynamics.

Currently working on the comparison of different time factors accompanying the between species, between patch and the combined response random numbers. To be updated ...

To check the preliminary results, please compile the linked executable:

*gfortran dvode_f90_m.f90 deterministic.f90 -o xx* (for deterministic trajectory calculations without Euler)

*gfortran det_euler.f90 -o xx* (for deterministic trajectory calculations with Euler)

*gfortran random.f90 xxxx.f90 -o xx* (for stochastic calculations - using Euler)

(please replace *gfortran* with your usual compiler id, *xxxx* by the file you would like to compile and let *xx* be the executable name.)

Steps to simulation:

1) Start with *deterministic.f90*. It generates system trajectories using *dvode_f90_m.f90*. It will give you a flavour of the species evolution.

2) Repeat the calculation with *det_euler.f90*. Compare the results as they should be identical.

3) ...
