# CFD
Contains cfd codes including my original Matlab codes from university and updated python codes.  
Codes:
  * rmannsol.m - Matlab Reimann-solver using linearized reimann equations, gives the exact solution at any spacetime point x/t  
  * reimann_solver.py - Python version, use the reimann_solver function to get values or use the output function to get a graph  
  * muscl.m - Matlab MUSCL-Hancock second-order solver, I remember it worked well when I last tried it but I have been having some issues with the Python port
  * blast.m - Matlab code which uses the MUSCL Hancock scheme to model a high-pressure/temperature region interacting with a low-temperaure/pressure region and a wall, i.e. the sequence of waves between a rapid combustion and a rigid wall.
