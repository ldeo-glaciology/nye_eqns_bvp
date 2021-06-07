%% Solving the Nye equations for the gradient in effective pressure and the gradient in discharge with Matlabs BVP code. 
% this is the first set of tests to check that the code is solving the
% equations 

% JK - April 17, 2020


%% setup x domain
xmesh = linspace(0,1,100);

%% define any channel size S
S = 0*xmesh+0.1;

%% define parameters
P.psi = 0.1;
P.e = 0.0034;
P.d = 1;
P.r = 0.9;
P.M = 0.001;

%% define boundary conditions on effective pressure at the lake Nl and the terminus Nt
NL = 0;  % effective pressure at the lake
Nt = 0;  % effective pressure at the terminus

%%  set some options
opts = bvpset('RelTol',0.00000001,'AbsTol',0.0000001,'Stats','on');

%% make an initial guess at the solutions
solinit= bvpinit(xmesh, @(x) Nyeguess(x,P));

%% compute the solutions: s5.y(1) = N;  s5.y(1) = Q; 
s5 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);

%% test 1a. check the derivatives
% check that the derivative computed by the function fed to the solver
% (Nye_NQ_1.m) with the solution as the input gives the same values as simply
% taking the gradient of the solution (i.e. are these values of Q and N
% actually solutions of these equations?).

N = s5.y(1,:);
Q = s5.y(2,:);
x = s5.x;
s = Nye_NQ(x,[N; Q],chan_size(xmesh,S,s5.x(1)),P);

dNdx = gradient(N)./gradient(x);
dQdx = gradient(Q)./gradient(x);

maxfracdiff(dNdx, s(1,:)) % maximum fractional difference is very small.
maxfracdiff(dQdx, s(2,:)) % maximum fractional difference is very small.

% (another sanity check: Nye_NQ_1 gives the same result as simply computing
% the righ-hand sides of the quesitons directly

S_x = chan_size(xmesh,S,x)';
rhs_1 = (abs(Q).*Q./S_x.^(8/3) - P.psi)/P.d;
rhs_2 = P.e*(P.r - 1)*abs(Q).^3./S_x.^(8/3) + P.e*S_x.*N.^3 + P.M;

maxfracdiff(dNdx, rhs_1) % maximum fractional difference is very small.

%% test 1b. does removing a term make a difference compared to making it eqaul zero?

P.e = 0;

s5_2 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);


s5_3 = bvp5c(@(x,y) Nye_NQ_test(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);

% answer: no
s5_2.stats.maxerr == s5_3.stats.maxerr  % it is exactly the same


%% test 2: For the BC at x=0 does it matter if you use N(0) = 0, or Q(0) computed by using N(0) = 0  ?
% reset the value of epsilon
P.e = 0.0034;

% Recompute the solution in the normal way with BCs on N defined 
s5_4 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,10,Nt),...
    solinit, opts);

x4 = s5_4.x;
N4 = s5_4.y(1,:);
Q4 = s5_4.y(2,:);

% define Q at the lake as the Q we got from above by prescribing NL
QL = Q4(1);   

% then solve the equations with QL defined 
s5_5 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_Qa_Nb(ya,yb,QL,Nt),...
    solinit, opts);

x5 = s5_5.x;
N5 = s5_5.y(1,:);
Q5 = s5_5.y(2,:);

maxfracdiff(N5,N4)  % the difference is small as long as epsilon is small

plot(x4,[N4;Q4],'b')
hold on
plot(x5,[N5;Q5],'g')


%% test 3: what about using a differnt S (but S remains uniform)?

S = 0.5+xmesh*0; 

s5_6 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,0,Nt),...
    solinit, opts);

x6 = s5_6.x;
N6 = s5_6.y(1,:);
Q6 = s5_6.y(2,:);

plot(x6,[N6;Q6],'r')    % it works and leads to a different solution 

%% test 3: What happens when S varies in space?
S = 0.5+0.45*sin(xmesh*pi*2); 
P.e = 0.0034;

P.psi=0
P.M = 0.1
s5_7 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,NL,Nt),...
    solinit, opts);

x7 = s5_7.x;
N7 = s5_7.y(1,:);
Q7 = s5_7.y(2,:);
hold on
plot(xmesh,S,'g'); hold on
plot(x7,[N7;Q7],'k')     % it effects the effective pressure and the discharge. In particular N is low around the constriction at x = 0.75. 


%% test 4: What happens when you use an initial guess very close (or equal) to the final result?


% re-compute a normal solution
S = 0.5+0*xmesh; 
P.e = 0.0034;

P.psi=0.1;
P.M = 0.001;
tic
solinit= bvpinit(xmesh, @(x) Nyeguess(x,P));
s5_8 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,0.1,0.11),...
    solinit, opts);
toc
x8 = s5_8.x;
N8 = s5_8.y(1,:);
Q8 = s5_8.y(2,:);
clf
plot(xmesh,S,'g'); hold on
plot(x8,[N8;Q8])

% there are two ways to feed the previous solution in as an initial guess
%  1 using a new guess function that takes as an input the previous
%  solution. 
tic
previous = s5_8;
solinit_test= bvpinit(xmesh, @(x) Nyeguess_previous(x,previous));
% this solves with a slightly different NL
s5_9 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,0.09,0.11),...
    solinit_test, opts);
toc
% it takes about 0.4s on my laptop, compared to about 1.3s using the standard initial
% guess form above:
tic
previous = s5_8;
% this solves with a slightly different NL
s5_9 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,0.09,0.11),...
    solinit, opts);
toc

%  2. feeding the bvpinit the old solution directly (this is simpler). 
tic
solinit_test= bvpinit(s5_8,[0 1]);
% this solves with a slightly different NL
s5_10 = bvp5c(@(x,y) Nye_NQ(x,y, chan_size(xmesh,S,x) ,P),...
    @(ya,yb) bc_N(ya,yb,0.09,0.11),...
    solinit_test, opts);
toc
% ^ this version is actually a little slower because it is using a higher
% resolution than the original xmesh that we prescribed at the top. 
