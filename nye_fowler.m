% Constants
lambda_l = 1;
Q_in = 1;
hL_pl1 = 1;
beta_r = 1;
delta = 1;
del_zeta = 1;

P.psi = 0.1;
P.e = 0.0034;
P.d = 1;
P.r = 0.9;
P.M = 0.001;

% Grid point sizing
n = 10; % space grid points
m = 10; % time grid points
del_s = 1; % space step size
del_t = 1; % time step size

% Initializing arrays to hold data
h = zeros(1,m); % Lake level is 1-D array evolving with time
% Other arrays evolve with time and along the channel
S = zeros(n,m);
Q = zeros(n,m);
N = zeros(n,m);

% Creating initial conditions
S(:,1) = 0.1*ones(n,1);

% Initializing boundary conditions
h(1,1) = 10;
NL = 0;  % effective pressure at the lake
Nt = 0;  % effective pressure at the terminus

% Looping to solve
for i = 1:del_t:m-1 % loop through time
    
    % Solve BVP for N_i, Q_i - adapted from Nye_BVP
    % Initialize array in space
    x_array = linspace(0,1,n);
    
    % Options for BVP solver - directly taken from Nye_BVP
    opts = bvpset('RelTol',0.00000001,'AbsTol',0.0000001,'Stats','on');
    
    % Initial Guess
    solinit= bvpinit(x_array, @(x) guess(x,P));
    
    %Solving
    s5 = bvp5c(@(x,y) Nye_NQ(x,y, x_array, S(:,i) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);
    N(:,i) = interp1(s5.x,s5.y(1,:),x_array);
    Q(:,i) = interp1(s5.x,s5.y(2,:),x_array);
    
    % Next use current N, Q, to get next step's S
    for j = 1:del_s:n % loop through space
        
        S(j,i+1) = S(j,i) + del_t*(abs(Q(j,i))^3/S(j,i)^(8/3) - S(j,i)*N(j,i)^3);
    end
    h(1,i+1) = h(1,i) + del_t*lambda_l*(Q_in-Q(1,i))/hL_pl1;
    
end

% Function for guessing
function y = guess(x,P)
y = [sin(x*pi)
    0+P.M*x];
end

% Function for N-Q equations
function dydx = Nye_NQ(x,y,xmesh,S_i,P)
S_x = interp1(xmesh,S_i,x);
dydx = [(y(2,:).*abs(y(2,:))./S_x.^(8/3) - P.psi)/P.d
    P.e*(P.r-1)*abs(y(2,:)).^3./S_x.^(8/3) + P.e*S_x.*y(1,:).^3 + P.M];
end

% Function for boundary conditions
function res = bc_N(ya,yb,NL,Nt) 
res = [ya(1)-NL
    yb(1)-Nt];
end
