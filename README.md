# Tester
(same as the wiki)

Simulation project for modeling metacommunity dynamics based on *Wang and Loreau* paper on "Biodiversity and ecosystem stability across scales in metacommunities" in *Ecology letters, Volume 19, Issue 5, May 2016, Pages 510–518*.

This compilaton/execution steps below are for Linux systems. For compiling and running these codes on Windows, please follow your Windows FORTRAN compiler procedure.

To check the preliminary results, please compile the linked executable:

(a) *gfortran dvode_f90_m.f90 deterministic.f90 -o xx* (for deterministic trajectory calculations without Euler)

(b) *gfortran det_euler.f90 -o xx* (for deterministic trajectory calculations with Euler)

(c) *gfortran random.f90 xxxx.f90 -o xx* (for stochastic calculations - using Euler)

(please replace *gfortran* with your usual compiler id, *xxxx* by the file you would like to compile and let *xx* be the executable name. Execution can be achieved by usual *./xx* command on the terminal.)

Steps to simulation:

(1) Start with *deterministic.f90*. Use procedure (a) above to generate the executable. It generates system trajectories using *dvode_f90_m.f90* in files *fort.11* and *fort.12*. It will give you a flavour of species evolution.

(2) Repeat the calculation with *det_euler.f90* (procedure (b)). The output is in files *fort.11* and *fort.12*. Compare the results as they should be identical to those obtained from step (1). Increase the evolution time and decrease the time step *h* if the results do not match. If issues still persist then please report, Thanks.
