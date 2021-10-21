% clear

% physical constants
g = 9.8;
L = 3.3e5;
rho_w =1000;
rho_i = 900;
K = 1e-24;
beta = 1000;

% geometry
H = 500;      % ice thickness
sigma = 0.01; % surface slope

% you have to choose a timescale
t0 = 60*60*100;   % 10 hour

% scales dependent on geometry
psi0 = rho_w*g*sigma;    % hydraulic gradient [Pa/m]   constant
x0 = 10000;    % length of channel [m]   constant
phi0 = 0.01;   % porosity of the glacier [dimensionless] a constant, but could be a variable in future versions
tau0 = rho_w*g*H*sigma;  % basal shear stress [Pa]

% scales dependent on climate
Q0 = 10;    % discharge [m^3/s], variable of the model

% scales that depend on other scales
M0 = Q0/x0;             % melt at surface and bed [m^2/s],input to the model
N0 = psi0*x0 ;           % effective pressue [Pa], variable of the model
m0 = Q0*psi0/L;          % melt rate along the channel walls [kg/s/m/s], variable of the model
S0 = m0/(rho_i*N0^3*K);  % channel cross-sectional area [m^2], variable of the model

% non-dimensional parameters
nu = S0*rho_i/t0/m0  ;             % channel area evolution constant
lambda = N0*phi0/t0/rho_w/g/M0 ;   % storage term evolution constant
delta = 1/phi0     ;               % porosity scale parameter
alpha = N0*x0/t0/tau0  ;           % porosity evolution parameter

ratio = rho_w*g*Q0/K/psi0^4/x0^5/phi0;   % nu/lamda first worked out analytically

nu/lambda;

N2head = @(N) (rho_i*g*H - N)/(rho_w*g);





