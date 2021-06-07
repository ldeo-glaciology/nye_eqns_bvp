function res = bc_Qa_Nb(ya,yb,QL,Nt) % boundary conditions
res = [ya(2)-QL
    yb(1)-Nt];
end