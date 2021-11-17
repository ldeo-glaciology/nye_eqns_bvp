% Function for solver
function [s_km,time_years,h_meters,Q_m3ps,S_m2,N_Pa,u_mpy] = nf_solver(grid_space,time_step)
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
s_to_y = 60*60*24*365.25;

% Scales
s0 = 10*10^3; % m
Q0 = 1500; % m^3/s
m0 = Q0*psi_0/L;
S0 = (f*rho_w*g*Q0^2/psi_0)^(3/8);
N0 = (K0*rho_i*S0*L/(psi_0*Q0))^(-1/3);
t0 = rho_i*S0*L/(psi_0*Q0); 
M0 = Q0/s0; 
h0 = 90;

% Model Params
P.e = s0*m0/(Q0*rho_i);
P.r = rho_i/rho_w;
P.d = N0/(s0*psi_0);

% Inputs
Q_in = 10/Q0;
hL_pl1 = 1;
beta_r = 0.9;
%u_raw = 100; % m/y
%u = u_raw*t0/(s_to_y*s0);
tau_raw = 100000; % Pa
tau = tau_raw/N0; 
C = 2*10^-20; % For sliding law
alpha = C*(N0^3)*t0/s0; % For sliding law

% Need to adjust
P.lambda = 3.2;
P.psi = 1;
P.M = 5*10^-4/M0; 

% Testing hack-y velocity max
umax = 1000/(s0*s_to_y/t0);
gamma = alpha*(tau^4)/umax;
% Grid point sizing
del_t = 0.1; % time step size

% Initializing arrays to hold data
h = zeros(1,time_step); % Lake level is 1-D array evolving with time
% Other arrays evolve with time and along the channel
S = zeros(grid_space,time_step);
Q = zeros(grid_space,time_step);
N = zeros(grid_space,time_step);
u = zeros(grid_space,time_step);
% Creating initial conditions
S(:,1) = 5*ones(grid_space,1)/S0;


% Initializing boundary conditions
h(1,1) = 1/3;
NL = beta_r*(1-h(1,1));  % effective pressure at the lake
Nt = rho_i*g*1/N0;  % effective pressure at the terminus - changed to constants
end_time = time_step-1;

% Initial step with guessing function
% Initialize array in space
x_array = linspace(0,1,grid_space);
P.psi_var = P.psi*(1-3*exp(-20.*x_array));

% Options for BVP solver - directly taken from Nye_BVP
opts = bvpset('RelTol',0.00000001,'AbsTol',0.0000001,'Stats','on');
solinit= bvpinit(x_array, @(x) Nyeguess(x,P));

% Solve with initial guess
s5 = bvp5c(@(x,y) Nye_NQ(x,y, x_array, S(:,1) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);
N(:,1) = interp1(s5.x,s5.y(1,:),x_array)';
Q(:,1) = interp1(s5.x,s5.y(2,:),x_array)';
u(:,1) = alpha*(tau^4)./(N(:,1)+gamma);
dSdx = gradient(S(:,1))./gradient(x_array)';
S(:,2) = S(:,1) + del_t.*(abs(Q(:,1)).^3./S(:,1).^(8/3) - S(:,1).*N(:,1).^3 - u(:,1).*dSdx);
h(1,2) = h(1,1) + del_t*P.lambda*(Q_in-Q(1,1))/hL_pl1;
NL = beta_r*(1-h(1,2));

% Looping to solve for rest of time
for i = 2:time_step-1 % loop through time

    % Guess
    solinit= bvpinit(x_array, @(x) [interp1(s5.x,s5.y(1,:),x)
                                        interp1(s5.x,s5.y(2,:),x)]);
    %Solving
    s5 = bvp5c(@(x,y) Nye_NQ(x,y, x_array, S(:,i) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);
    N(:,i) = interp1(s5.x,s5.y(1,:),x_array)';
    Q(:,i) = interp1(s5.x,s5.y(2,:),x_array)';
    u(:,i) = alpha*(tau^4)./(N(:,i)+gamma);
    % Next use current N, Q, S to get next step's S and h
    dSdx = gradient(S(:,i))./gradient(x_array)';
    S(:,i+1) = S(:,i) + del_t.*(abs(Q(:,i)).^3./S(:,i).^(8/3) - S(:,i).*N(:,i).^3 - u(:,i).*dSdx);
    h(1,i+1) = h(1,i) + del_t*P.lambda*(Q_in-Q(1,i))/hL_pl1;
    if h(1,i+1) <= 0 || S(1,i+1) <= 0 % Adding error handling for if channel closes/lake drains
        end_time = i+1;
        break
    end
    NL = beta_r*(1-h(1,i+1));
    
end

% Convert to scale properly
s_km = s0*x_array/1000;
time_years = t0*del_t*[1:end_time]/s_to_y;
h_meters = h0*h(1,1:end_time);
Q_m3ps = Q0*Q(:,1:end_time);
S_m2 = S0*S(:,1:end_time);
N_Pa = N0*N(:,1:end_time);
u_mpy = u*s0*s_to_y/t0;
end


