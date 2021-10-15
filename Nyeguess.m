% Function for initial guess for N (y(1)) and Q (y(2))
function y = Nyeguess(x,P)
y = [sin(x*pi)
    0.001+(P.M)*x]; % try making first 0 nonzero, but small
end
