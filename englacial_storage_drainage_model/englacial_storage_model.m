
scales

%% space domain
dx = 0.01;
X=1;
x = 0:dx:X;
Lx = length(x);

%% time domain
dt = 0.000001;
T=10;
t = 0:dt:T;
Lt = length(t);

%% initial conditions

S = x*0 + 1;
N = (X-x)/X*0.1;

%% psi
psi = 1;

%% M
M = 0.1;

%% phi
phi = 1;

Qtop = 0;
Qbottom = 0.1;

for ii = 2:Lt
    
    N_old = N;
    S_old = S;
    if any(imag(S)>0)
        pause(0.001)
    end
    S_stag = (S_old(2:end) + S_old(1:end-1))/2;
    PSI = (psi + (N_old(2:end) - N_old(1:end-1))/dx);
    Q_stag = sign(PSI).*sqrt(S_stag.^(8/3).*abs(PSI));  % Q on the staggered grid
    
    N(2:end-1) = N_old(2:end-1) + dt/lambda*( (Q_stag(2:end) - Q_stag(1:end-1))/dx  -  M)/phi;
    N(1) = 0; %N(2) + dx*(-psi + Qtop^2/S(1)^(8/3));
    N(end) = 0.1; %N(end-1) + dx*(psi - Qbottom^2/S(1)^(8/3));
    
    Q = 0*S;
    PSI = psi + (N(3:end) - N(1:end-2))/(2*dx);
    if any(imag(PSI)>0)
        pause(0.001)
    end
    Q(2:end-1) = sign(PSI).*sqrt(S(2:end-1).^(8/3).*(abs(PSI)));
    Q(1) = Qtop; % consistent with the N gradient imposed above.
    Q(end) = Qbottom;
    S = S_old + dt/nu*(abs(Q).^3 - S.*N.^3);
    
    
    if rem(ii,10000)==0
        plot(x,S) %x,N,x,Q,
        %legend('N','Q','S')
        %ylim(yl)
        drawnow
        %     pause
         t(ii)
    end
end
