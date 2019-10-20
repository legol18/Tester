# Tester
(same as the wiki)

Simulation project for modeling metacommunity dynamics based on *Wang and Loreau* paper on "Biodiversity and ecosystem stability across scales in metacommunities" in *Ecology letters, Volume 19, Issue 5, May 2016, Pages 510–518*.

This compilaton/execution steps below are for Linux systems. For compiling and running these codes on Windows, please follow your Windows FORTRAN compiler procedure

To check the preliminary results, please compile the linked executable:

(a) *gfortran dvode_f90_m.f90 deterministic.f90 -o xx* (for deterministic trajectory calculations without Euler)

(b) *gfortran det_euler.f90 -o xx* (for deterministic trajectory calculations with Euler)

(c) *gfortran random.f90 xxxx.f90 -o xx* (for stochastic calculations - using Euler)

(please replace *gfortran* with your usual compiler id, *xxxx* by the file you would like to compile and let *xx* be the executable name. Execution can be achieved by usual *./xx* command on the terminal.)

Steps to simulation:

You can play around with all the codes and procedures as you like. 

To generate the diversity vs variability plot from the paper, you need the following routines in the same folder, namely:

(a) latest.f90, 

(b) random.f90, 

(c) plotting.gpl (gnuplot file, needs gnuplot installed to work), 

(d) exec.sh (Primary execution file. Please enable appropriate permissions using chmod before executing).

Please report issues, Thanks.
