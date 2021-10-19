function dydx = Nye_NQ(x,y,xmesh,S_i,P) % equation to solve
% encodes the dN/dx and dQ/dx equations from page 28 of Kingslake (2013) in
% the format that bvp5c needs. 

% J. Kingslake, April 2020. 
% Updated to fit George Lu's N-F solver, October 2021

S_x = interp1(xmesh,S_i,x);
psi_x = interp1(xmesh,P.psi_var,x);

S83 = S_x.^(8/3);
N = y(1,:);
Q = y(2,:);
dydx = [(Q.*abs(Q)./S83 - psi_x )/P.d % change between constant/variable psi
    P.e*(P.r-1)*abs(Q).^3./S83 + P.e*S_x.*N.^3 + P.M]; % Use just P.M for efficiency
end
