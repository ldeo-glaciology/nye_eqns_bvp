# nye_eqns_bvp
Solving the Nye-fowler equations describing subglacial channelized flow. 

Nye_BVP.m performs tests using Matlab's built-in boundary value to solve coupled for effective pressure, N, and discharge, Q.

All other functions are called by Nye_BVP.m. For example, Nye_NQ.m encodes the PDEs in the format required by matlab's BVP solvers. 
