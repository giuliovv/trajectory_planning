clear all
close all
clc

%% Parameters of car and circuit
w=1.992;
e=0.40;
approx_circuit_width = 6;
rel_e = e/approx_circuit_width;

%initialization of kr and kl, tuning of their ratio
kl= 0;
kr=1e7;

%gamma in case of border constraint
gamma = 0;

%track borders
[rightpoints, leftpoints] = read_track("../tracks/Monza.csv");
n = height(rightpoints);

%% Initialization of the coefficient alpha
alpha0 = 0.5*ones(n, 1);

%% Linear equality constraint parameters
A               =   [];
b               =   [];

%% Linear inequality constraint parameters
C               =   [eye(n); -eye(n)];
d               =   [zeros(n,1)+rel_e; -(ones(n,1)-rel_e)];

%% Optimization parameters
myoptions   =   myoptimset;
myoptions.ls_beta       = 0.3;        
myoptions.ls_c          = 0.1;
myoptions.gradmethod    = 'UP';
myoptions.graddx        = eps^(1/3);
myoptions.nitermax      = 1e2;
myoptions.tolfun        =	1e-12;
myoptions.Hessmethod    = 'GN';
myoptions.GN_funF       = @(alpha)mycostfunction_constrained_GN(n,rightpoints,leftpoints,alpha, kl, kr, gamma);
myoptions.GN_sigma      =	1e2;

myoptions.ls_nitermax   =	30;

%% Optimization routine
%different coordinates to pass through
%cost function

[Ustar,fxstar,niter,exitflag,xsequence] =myfmincon(@(x)mycostfunction_constrained(n,rightpoints,leftpoints,alpha, kl, kr, gamma),alpha0,A,b,C,d,0,1,myoptions);

%% Plot
 [x,y] = mycoordinates(n,rightpoints,leftpoints,Ustar);
hold on
plot(leftpoints(:, 1), leftpoints(:, 2), 'k')
plot(rightpoints(:, 1), rightpoints(:,2), 'k')
plot(x, y, 'r')
hold off

%% Laptime
time = laptime(x, y, n, 100/2.9, 11)