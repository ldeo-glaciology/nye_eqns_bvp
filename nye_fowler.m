% Constants
rho_w = 1000; % kg/m^3 (water density)
rho_i = 900; %kg/m^3 (ice density)
g = 10; % m/s^2 (gravitational acceleration)
L = 3.3*10^5; % J/kg
K0 = 10^-24; % Pa^-3 s^-2
phi_s = 0.01; % surface slope
psi_0 = rho_w*g*phi_s;
n_prime = 0.1; % m^-1/3 s (hydraulic roughness)
f = 6.6*n_prime^2; 

% Scales
s0 = 10*10^3; % m
Q0 = 1500; % m^3/s
m0 = Q0*psi_0/L;
S0 = (f*rho_w*g*Q0^2/psi_0)^(3/8);
N0 = (K0*rho_i*S0*L/(psi_0*Q0))^(-1/3);
t0 = rho_i*S0*L/(psi_0*Q0); 
M0 = Q0/s0; 

% Model Params
e = s0*m0/(Q0*rho_i);
r = rho_i/rho_w;
d = N0/(s0*psi_0);
% to be continued

% Inputs
Q_in = 10/Q0;
hL_pl1 = 1;
beta_r = 0.9;

P.lambda = 3.2;
P.psi = 1;
P.e = 0.0034;
P.d = 1;
P.r = 0.9;
P.M = 0.00; % I can't change to zero without error - might have to do w guessing

% Grid point sizing
n = 10; % space grid points
m = 5000; % time grid points
del_s = 0.01; % space step size
del_t = 0.1; % time step size

% Initializing arrays to hold data
h = zeros(1,m); % Lake level is 1-D array evolving with time
% Other arrays evolve with time and along the channel
S = zeros(n,m);
Q = zeros(n,m);
N = zeros(n,m);

% Creating initial conditions
S(:,1) = 0.01*ones(n,1);

P.psi_var = P.psi*(1-3*exp(-20*S(:,1)./S(:,1)));

% Initializing boundary conditions
h(1,1) = 1/3;
NL = beta_r*(1-h(1,1));  % effective pressure at the lake
Nt = 0;  % effective pressure at the terminus
end_time = m-1;
% Looping to solve
for i = 1:m-1 % loop through time
    
    % Solve BVP for N_i, Q_i - adapted from Nye_BVP
    % Initialize array in space
    x_array = linspace(0,1,n);
    
    % Options for BVP solver - directly taken from Nye_BVP
    opts = bvpset('RelTol',0.00000001,'AbsTol',0.0000001,'Stats','on');
    
    % Initial Guess
    if i == 1
        solinit= bvpinit(x_array, @(x) guess(x,P));
    else % Use previous solution to make guess
        solinit= bvpinit(x_array, @(x) [interp1(s5.x,s5.y(1,:),x)
                                        interp1(s5.x,s5.y(2,:),x)]);
    end
    
    %Solving
    s5 = bvp5c(@(x,y) Nye_NQ(x,y, x_array, S(:,i) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);
    N(:,i) = interp1(s5.x,s5.y(1,:),x_array);
    Q(:,i) = interp1(s5.x,s5.y(2,:),x_array);
    
    % Next use current N, Q, S to get next step's S and h
        
    S(:,i+1) = S(:,i) + del_t.*(abs(Q(:,i)).^3./S(:,i).^(8/3) - S(:,i).*N(:,i).^3);

    h(1,i+1) = h(1,i) + del_t*P.lambda*(Q_in-Q(1,i))/hL_pl1;
    if h(1,i+1) <= 0 || S(1,i+1) <= 0 % Adding error handling for if channel closes/lake drains
        end_time = i+1;
        break
    end
    NL = beta_r*(1-h(1,i+1));
    P.psi_var = P.psi*(1-3*exp(-20*S(:,i)./S(:,1)));
    
end
%% Plotting

% Making plots
s_to_y = 60*60*24*365.25;
figure('Name','Evolution of Drainage');
subplot(3,1,1);
plot(t0*del_t*[1:end_time]/s_to_y,100*h(1,1:end_time));
xlabel('Time (years)');
ylabel('Height');
title('Lake height over time');

subplot(3,1,2);
plot(t0*del_t*[1:end_time]/s_to_y,Q0*Q(1,1:end_time),'DisplayName','Lake Exit');
hold on;
plot(t0*del_t*[1:end_time]/s_to_y,Q0*Q(n,1:end_time),'DisplayName','Channel Exit');
xlabel('Time (years)');
ylabel('Flow Rate (m^3 s^-1');
title('Flow over time');
legend('Location','northwest');

subplot(3,1,3);
plot(t0*del_t*[1:end_time]/s_to_y,S0*S(1,1:end_time),'DisplayName','Lake Exit'); % now plotting area at lake
hold on;
plot(t0*del_t*[1:end_time]/s_to_y,S0*S(n,1:end_time),'DisplayName','Channel Exit');
xlabel('Time (years)');
ylabel('Channel Area (m^2)');
title('Channel area over time');
legend('Location','northwest');

%% Functions
% Function for guessing
function y = guess(x,P)
y = [sin(x*pi)
    0.001+(P.M)*x]; % try making first 0 nonzero, but small
end

% Function for N-Q equations
function dydx = Nye_NQ(x,y,xmesh,S_i,P)
S_x = interp1(xmesh,S_i,x);
psi_x = interp1(xmesh,P.psi_var,x);
dydx = [(y(2,:).*abs(y(2,:))./S_x.^(8/3) - P.psi )/P.d
    P.e*(P.r-1)*abs(y(2,:)).^3./S_x.^(8/3) + P.e*S_x.*y(1,:).^3 + P.M];
end

% Function for boundary conditions
function res = bc_N(ya,yb,NL,Nt) 
res = [ya(1)-NL
    yb(1)-Nt];
end


