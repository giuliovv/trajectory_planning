clear all
close all
clc

%% Parameters of car and circuit
w=1.992;
e=0.40;

%track borders
[rightpoints, leftpoints] = read_track("tracks/Monza.csv");
n = height(rightpoints);

%% Initialization of the coefficient alpha
alpha0 = 0.5*zeros(n);

%% Optimization parameters
myoptions   =   myoptimset;
myoptions.ls_beta       = 0.3;        
myoptions.ls_c          = 0.1;
myoptions.gradmethod    = 'IM';
myoptions.graddx        = eps^(1/3);
myoptions.nitermax      = 5e2;
myoptions.Hessmethod    = 'SD';
%myoptions.GN_funF       = TODO

%% Optimization routine
%different coordinates to pass through
%cost function
J = mycostfunction(n,rightpoints,leftpoints,alpha0);
J

[Ustar,fxstar,k,exitflag,xsequence] = myfminunc(@(alpha)mycostfunction(n,rightpoints,leftpoints,alpha),alpha0,myoptions);