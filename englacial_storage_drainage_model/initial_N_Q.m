n=1
tol = 1e-7;
N_old = -9999+N;
dt_temp = 3.5/t0; 
while max(abs(N_old - N)) > tol
    N_old = N;
    
    N(2:end-1) = N(2:end-1) + dt_temp/lambda*( (Q(2:end) - Q(1:end-1))/dx  -  M)/phi;  % on the main grid
    
    
    PSI = (psi + (N(2:end) - N(1:end-1))/dx);   % PSI on the staggered grid
    
    Q = sign(PSI).*sqrt(S.^(8/3).*abs(PSI));
    
%     if rem(n,5000)==0
%         if ~ishandle(1); figure(1); end
%         set(0, 'CurrentFigure', 1)
%         plot(x,N,x_stag,Q,x_stag,S) ;  legend('N','Q','S')
%         drawnow
%         n
%         
%     end
%     n=n+1;
    
end