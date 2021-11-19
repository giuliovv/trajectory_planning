  function J= mycostfunction(n,rightpoints,leftpoints,alpha)

% MYCOSTFUNCTION computes the value of the cost function. 
%  INPUTS:  n             = number of segments dividing the track
%           rightpoints   = coordinates of the inner points of the track
%           leftpoints    = coordinates of the outer points of the track
%           alpha         = optimization variable, needed to express where
%                          the car has to stay inside the track
%
%  OUTPUTS: J             = value of the cost function

%initialization of kr and kl, tuning of their ratio
kl=1;
kr=1;


[x,y] = mycoordinates(n,rightpoints,leftpoints,alpha);

%cost function
Jl=0;
for idx = 2:n-2
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    Jl=Jl+ds^2;
end
Jr=0;
for idx = 3:n-2
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    dtheta=atan((y(idx+1)-y(idx))/(x(idx+1)-x(idx)))-atan((y(idx)-y(idx-1))/(x(idx)-x(idx-1)));
    rho=dtheta/ds;
    Jr=Jr+rho^2;
end
J=kl*Jl+kr*Jr;

