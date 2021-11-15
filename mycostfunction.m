function J= mycostfunction(n,x,y)

% MYCOSTFUNCTION computes the value of the cost function. 
%  INPUTS:  n             = number of segments dividing the track
%           x             = values of x at each step
%           y             = values of y at each step
%
%  OUTPUTS: J             = value of the cost function

%initialization of kr and kl, tuning of their ratio
kl=1;
kr=1;

%cost function
Jl=0;
for idx = 0:n-1
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    Jl=Jl+ds^2;
end
Jr=0;
for idx = 1:n-1
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    dtheta=atan((y(idx+1)-y(idx))/(x(idx+1)-x(idx)))-atan((y(idx)-y(idx-1))/(x(idx)-x(idx-1)));
    rho=dtheta/ds;
    Jr=Jr+rho^2;
end
J=kl*Jl+kr*Jr;

