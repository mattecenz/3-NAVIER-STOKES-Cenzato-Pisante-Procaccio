# 3-NAVIER-STOKES-Cenzato-Pisante-Procaccio
Project of the course Numerical Methods for Partial Differential Equations of Politecnico di Milano.


### GENERAL INFORMATIONS
This project aims to solve the unsteady, incompressible Navier-Stokes equations using the finite element method. The focus is on simulating the benchmark problem "flow past a cylinder" in two or three dimensions.

### TO RUN IT
To compile and run the project these are the steps that need to be followed:

+ if not already there create a folder build
+ move inside the folder
+ load the dealii modules with this command:`module load gcc-glibc dealii`
+ execute: `cmake ..`
+ execute `make`
+ if you want to run the 2D test execute `./navier_stokes2D`
+ if you want to run the 3D test execute `./navier_stokes3D`

Both these tests can be run also in parallel with MPI.
