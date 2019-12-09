# pipe_channel_flow repository
Fortran code for Direct Numerical Simulation of pipe and channel flows. Code implements a finite-difference method for incompressible Navier–Stokes equations solving described in https://doi.org/10.1016/j.jcp.2006.01.036. 

How is this repository?
=======================
Primary branches: 
+ pipe - simplest version of code for pipe flow time integration and support tools.
+ pipeSym - code for pipe flow time integration with symmetry reduction (reflection and rotation symmetry) and support tools.
+ pipeMPI - MPI version of code from pipe branch. 
+ pipeSymMPI - MPI version of code from pipeSym branch.
+ ductSym - code for channel flow simulation with additional symmetries.
+ ductSymMPI - MPI version of code from ductSym branch.
+ pytools - python scripts for experiments processing
+ master - not in use

Other branches with names like (primary_branch_name)-(feature) produced from primary branches and contain special cases of code.
