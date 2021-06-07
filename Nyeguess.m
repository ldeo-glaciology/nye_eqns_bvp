function g = Nyeguess(x,P) % initial guess for N (y(1)) and Q (y(2))
g = [sin(x*pi)
    0+P.M*x];
end