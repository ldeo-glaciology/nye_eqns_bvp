function dydx = Nye_NQ(x,y,S_x,P) % equation to solve
% encodes the dN/dx and dQ/dx equations from page 28 of Kingslake (2013) in
% the format that bvp5c needs. 

% J. Kingslake, April 2020. 

dydx = [(y(2,:).*abs(y(2,:))./S_x.^(8/3) - P.psi)/P.d
    P.e*(P.r-1)*abs(y(2,:)).^3./S_x.^(8/3) + P.e*S_x.*y(1,:).^3 + P.M];
end