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

(3) (issues) Now after checking the calculations for the deterministic cases, compile the file *astoc_sim1.f90* as per procedure (c) mentioned above and execute. The program will ask for an input which corresponds to the "between patch correlation" value. Based on the paper, pick a value in *(-0.8,0.8)* and check the output trajectories stored in files *fort.11* and *fort.12*. Please note that the program will also print two rows as screen outputs. Here the first row corresponds to the state of the system (equilibrium in this system's case) after removing the transients and the second row corresponds to one possible state of the system after the noisy evolution. I placed these to check for numerical issues during the evolution. Corresponding lines can be removed without changing the required output. 

(4) (issues) Copy *astoc_sim2.f90*, *random.f90* and *exec.sh* is a new folder, give *exec.sh* execution permission and run the file. Data files with different biodiversity parameters for different "between patch correlation" values are generated, namely for *-0.8*, *0.0* and *0.8* ...
