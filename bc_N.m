% Function for boundary conditions
function res = bc_N(ya,yb,NL,Nt) 
res = [ya(1)-NL
    yb(1)-Nt];
end