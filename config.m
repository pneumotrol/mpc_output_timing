close all;
clear;
clc;

%% simulatoin definition
% dt = 0.02; % stable
dt = 0.2; % unstable
tf = 4;
x0 = [1;0];

%% system definition (mass-damper-spring)
m = 1;
d = 1;
k = 1;

A = [
    0,1;
    -k/m,-d/m;
    ];

B = [
    0;
    1/m;
    ];

C = [
    1,0;
    ];

D = 0;

% dimensions
l = size(C,1); % output
m = size(B,2); % input
n = size(A,2); % state

% continuous state space
sysc= ss(A,B,C,D);

% discrete state space
sysd = c2d(sysc,dt,"zoh");

%% conventional MPC
N = ceil(tf/dt);
Q = diag([1]);
R = diag(0.01);
c = mpc_coeffs(sysd,struct("N",N,"Q",Q,"R",R));
MPC.K = (c.G'*c.Q*c.G + c.R)\(c.G'*c.Q);
MPC.F = c.F;
MPC.I = repmat(eye(l),N,1);

%% MPC w/ delay compensation
N = ceil(tf/dt);
Q = diag([1]);
R = diag(0.01);
c = mpc_coeffs(sysd,struct("N",N,"Q",Q,"R",R));
F = [c.F,c.G(:,1:m)];
F = F(l+1:end,:);
G = c.G(l+1:end,m+1:end);
Q = c.Q(l+1:end,l+1:end);
R = c.R(m+1:end,m+1:end);

MPC_compensated.K = (G'*Q*G + R)\(G'*Q);
MPC_compensated.F = F;
MPC_compensated.I = repmat(eye(l),N-1,1);

%% simulation
sim("model.slx");
