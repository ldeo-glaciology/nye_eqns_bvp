function g = Nyeguess_previous(x,previous) % initial guess for N (y(1)) and Q (y(2))
g = nan(1,length(x));



g(1,:) = interp1(previous.x,previous.y(1,:),x)+0*rand(1,length(x));
g(2,:) = interp1(previous.x,previous.y(2,:),x)+0*rand(1,length(x));

end