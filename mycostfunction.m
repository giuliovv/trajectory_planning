  function [j, gradient]= mycostfunction(n,rightpoints,leftpoints,alpha)

% MYCOSTFUNCTION computes the value of the cost function. 
%  INPUTS:  n             = number of segments dividing the track
%           rightpoints   = coordinates of the inner points of the track
%           leftpoints    = coordinates of the outer points of the track
%           alpha         = optimization variable, needed to express where
%                          the car has to stay inside the track
%
%  OUTPUTS: J             = value of the cost function
%               gradient     = value of the cost function gradient

%initialization of kr and kl, tuning of their ratio
kl=1;
kr=1;


[x,y] = mycoordinates(n,rightpoints,leftpoints,alpha);

%cost function
Jl=0;
jld = 0;
for idx = 2:n-2
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    Jl=Jl+ds^2;
    dsd = (alpha(idx)*(leftpoints(idx+1,1)-rightpoints(idx+1,1)-leftpoints(idx,1)+rightpoints(idx,1))+rightpoints(idx+1,1)-rightpoints(idx,1))* 2*(leftpoints(idx+1,1)-leftpoints(idx,1)-rightpoints(idx+1,1)+rightpoints(idx,1))+(alpha(idx)*(leftpoints(idx+1,2)-rightpoints(idx+1,2)-leftpoints(idx,2)+rightpoints(idx,2))+rightpoints(idx+1,2)-rightpoints(idx,2))*2*(leftpoints(idx+1,2)-leftpoints(idx,2)-rightpoints(idx+1,2)+rightpoints(idx,2));
    Jld=Jld+dsd;
end
Jr=0;
for idx = 3:n-2
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    dtheta=atan((y(idx+1)-y(idx))/(x(idx+1)-x(idx)))-atan((y(idx)-y(idx-1))/(x(idx)-x(idx-1)));
    rho=dtheta/ds;
    Jr=Jr+rho^2;
    dsd =(-((-leftpoints(idx-1, 2) + leftpoints(idx, 2) + rightpoints(idx-1, 2) - rightpoints(idx, 2))/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1)) + (-alpha(idx)*(leftpoints(idx-1, 2) - rightpoints(idx-1, 2)) + alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) - rightpoints(idx-1, 2) + rightpoints(idx, 2))*(leftpoints(idx-1, 1) - leftpoints(idx, 1) - rightpoints(idx-1, 1) + rightpoints(idx, 1))/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1))^2)/(1 + (-alpha(idx)*(leftpoints(idx-1, 2) - rightpoints(idx-1, 2)) + alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) - rightpoints(idx-1, 2) + rightpoints(idx, 2))^2/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1))^2) + ((leftpoints(idx+1, 2) - leftpoints(idx, 2) - rightpoints(idx+1, 2) + rightpoints(idx, 2))/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1)) + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))*(-leftpoints(idx+1, 1) + leftpoints(idx, 1) + rightpoints(idx+1, 1) - rightpoints(idx, 1))/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2)/(1 + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))^2/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2))/sqrt((alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2 + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))^2) + (-1.0/2.0*(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))*(2*leftpoints(idx+1, 1) - 2*leftpoints(idx, 1) - 2*rightpoints(idx+1, 1) + 2*rightpoints(idx, 1)) - 1.0/2.0*(alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))*(2*leftpoints(idx+1, 2) - 2*leftpoints(idx, 2) - 2*rightpoints(idx+1, 2) + 2*rightpoints(idx, 2)))*(atan((alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))) - atan((-alpha(idx)*(leftpoints(idx-1, 2) - rightpoints(idx-1, 2)) + alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) - rightpoints(idx-1, 2) + rightpoints(idx, 2))/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1))))/((alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2 + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*((leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))^2))^(3.0/2.0);
    jrd = jrd+dsd;
end
J=kl*Jl+kr*Jr;
gradient = kl*Jld+kr*Jrd;

