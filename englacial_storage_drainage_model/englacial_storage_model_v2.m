clear
%close all
scales
tic

%% space domain
dx = 0.05;
X=1;
x = 0:dx:X;
x_stag = (x(2:end)+x(1:end-1))/2;
Lx = length(x);
HalfLx = round(Lx/2);

%% time domain
dt = 0.0002;
T=10000;
t = 0:dt:T;
Lt = length(t);
HalfLt = round(Lt/2);
%% psi
psi = 1;

%% S
S = x_stag*0 + 0.1;

%% M
M = 0.1;

%% BCs

Ntop = 0.3;   %+t*0;%+0.05*sin(2*pi*t*t0/2000);
Nbottom = 0;

%% phi    porosity (nondimensional)
phi = 50;


%% initial condition
% a random profile that obeys the boundary conditions. 
N = Ntop(1) + x*(Nbottom-Ntop(1))/X; %+ 0.2*sin(x*pi);
PSI = (psi + (N(2:end) - N(1:end-1))/dx);   % PSI on the normal grid, but not the 1st and last point
Q = sign(PSI).*sqrt(S.^(8/3).*abs(PSI)) ;
initial_N_Q

%% pre-allocate output arrays
Qout = t*nan;
Sout = t*nan;
Nout = t*nan;

Qout(1) = Q(HalfLx);
Sout(1) = S(HalfLx);
Nout(1) = N(HalfLx);


for ii = 2:Lt
    N_old = N;
    Q_old = Q;
    % update N
    N(2:end-1) = N(2:Lx-1) + dt/lambda * ( (Q(2:(Lx-1)) - Q(1:(Lx-2)))/dx  -  M)/phi;  % on the main grid
%    N(1) = Ntop;
    
    N_stag = (N(2:Lx) + N(1:(Lx-1)))/2;
    
    % update S
    S = S + dt/nu*(abs(Q_old).^3 - S.*N_stag.^3);
    
    PSI = (psi + (N(2:Lx) - N(1:(Lx-1)))/dx);   % PSI on the staggered grid
    
    Q = sign(PSI).*sqrt(S.^(8/3).*abs(PSI));  % G on the staggered grid
    
    
     if rem(ii,50000)==0
%         if ~ishandle(1); figure(1); end
%         set(0, 'CurrentFigure', 1)
%         plot(x,N,x_stag,Q,x_stag,S) ;  legend('N','Q','S')
%         %axis([0 1 0 1])
%         %         if ~ishandle(2); figure(2); end
%         %         set(0, 'CurrentFigure', 2)
%         %         plot(x,(N-N_old)/dt) ; legend('dN/dt')
%         
%         %         if ~ishandle(3); figure(3); end
%         %         set(0, 'CurrentFigure', 3)
%         %         hold on
%         %         plot(t(ii)*t0,N0/t0/rho_w/g*(N(round(Lx/2))-N_old(round(Lx/2)))/dt,'g*') ; legend('dN(center,t)/dt')
%         if ~ishandle(4); figure(4); end
%         set(0, 'CurrentFigure', 4)
%         plot(N(4),Q(4),'.') ;  title('N vs Q'); hold on
%         %         axis([0 1 0 1])
%         
%         
%         %ylim(yl)
%         drawnow
%         %pause
         disp([num2str(t(ii)*t0/60/60/24) ' out of ' num2str(T*t0/60/60/24) ' days'])
     end

    Qout(ii) = Q(HalfLx);
    Sout(ii) = S(HalfLx);
    Nout(ii) = N(HalfLx);

end
toc

hold off
plot(t,Nout)
