clear all
clean all
clc

%parameters of car and circuit
w=1.992;
e=0.40;

%qui dobbiamo introdurre il circuito (abbiamo bisogno dei parametri n e
%destra e sinistra), basta mettere il comando che fa partire il circuito?

%initialization of the coefficient alpha
alpha = 0.5*zeros(n);

%different coordinates to pass through
[x,y] = mycoordinates(n,rightpoints,leftpoints,alpha);

%cost function
J = mycostfunction(n,x,y);