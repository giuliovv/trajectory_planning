clear all
close all
clc

%parameters of car and circuit
w=1.992;
e=0.40;

%track borders
[rightpoints, leftpoints] = read_track("tracks/Monza.csv");
n = height(rightpoints);

%initialization of the coefficient alpha
alpha = 0.5*zeros(n);

%different coordinates to pass through
[x,y] = mycoordinates(n,rightpoints,leftpoints,alpha);

%cost function
J = mycostfunction(n,x,y);