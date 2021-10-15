function dydx = Nye_NQ(x,y,xmesh,S_i,P) % equation to solve
% encodes the dN/dx and dQ/dx equations from page 28 of Kingslake (2013) in
% the format that bvp5c needs. 

% J. Kingslake, April 2020. 
% Updated to fit George Lu's N-F solver, October 2021

S_x = interp1(xmesh,S_i,x);
psi_x = interp1(xmesh,P.psi_var,x);
dydx = [(y(2,:).*abs(y(2,:))./S_x.^(8/3) - psi_x )/P.d % change between constant/variable psi
    P.e*(P.r-1)*abs(y(2,:)).^3./S_x.^(8/3) + P.e*S_x.*y(1,:).^3 + P.M];
end
