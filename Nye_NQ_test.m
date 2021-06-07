function dydx = Nye_NQ_test(x,y,S_x,P) % equation to solve
% encodes the dN/dx and dQ/dx equations from page 28 of Kingslake (2013) in
% the format that bvp5c needs. 

% test: removing the epsilon terms 

% J. Kingslake, April 2020. 

dydx = [(y(2,:).*abs(y(2,:))./S_x.^(8/3) - P.psi)/P.d
        P.M];
end