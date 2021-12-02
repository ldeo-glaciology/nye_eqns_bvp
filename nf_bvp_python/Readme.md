# Nye Fowler Solver
This set of code solves modified Nye-Fowler equations to simulate the drainage of a glacier lake, and its (hopefully soon-to-be) coupled glacier dynamics. The details on the equations and each version can be found at: https://www.overleaf.com/read/hzvvwvjkfwfj

A condensed version history for nf_bvp.ipnb can be found below:

- v1.0: Python Version of N-F solver solving the base equations
- v2.1: Added a velocity adjustment using a basic sliding law 
- v3.1.1: Added a psi defined by ice thickness (with heavy constraints)