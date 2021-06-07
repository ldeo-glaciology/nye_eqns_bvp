function res = bc_N(ya,yb,NL,Nt) % boundary conditions
res = [ya(1)-NL
    yb(1)-Nt];
end