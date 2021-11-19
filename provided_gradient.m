  function J= provided_gradient(n,rightpoints,leftpoints,alpha)

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

%cost function
Jld=0;
for idx = 2:n-2
     dsd = (alpha(idx)*(leftpoints(idx+1,1)-rightpoints(idx+1,1)-leftpoints(idx,1)+rightpoints(idx,1))+rightpoints(idx+1,1)-rightpoints(idx,1))* 2*(leftpoints(idx+1,1)-leftpoints(idx,1)-rightpoints(idx+1,1)+rightpoints(idx,1))+(alpha(idx)*(leftpoints(idx+1,2)-rightpoints(idx+1,2)-leftpoints(idx,2)+rightpoints(idx,2))+rightpoints(idx+1,2)-rightpoints(idx,2))*2*(leftpoints(idx+1,2)-leftpoints(idx,2)-rightpoints(idx+1,2)+rightpoints(idx,2));
    Jld=Jld+dsd;
end
Jr=0;
for idx = 3:n-2
   dsd =(-((-leftpoints(idx-1, 2) + leftpoints(idx, 2) + rightpoints(idx-1, 2) - rightpoints(idx, 2))/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1)) + (-alpha(idx)*(leftpoints(idx-1, 2) - rightpoints(idx-1, 2)) + alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) - rightpoints(idx-1, 2) + rightpoints(idx, 2))*(leftpoints(idx-1, 1) - leftpoints(idx, 1) - rightpoints(idx-1, 1) + rightpoints(idx, 1))/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1))^2)/(1 + (-alpha(idx)*(leftpoints(idx-1, 2) - rightpoints(idx-1, 2)) + alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) - rightpoints(idx-1, 2) + rightpoints(idx, 2))^2/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1))^2) + ((leftpoints(idx+1, 2) - leftpoints(idx, 2) - rightpoints(idx+1, 2) + rightpoints(idx, 2))/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1)) + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))*(-leftpoints(idx+1, 1) + leftpoints(idx, 1) + rightpoints(idx+1, 1) - rightpoints(idx, 1))/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2)/(1 + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))^2/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2))/sqrt((alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2 + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))^2) + (-1.0/2.0*(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))*(2*leftpoints(idx+1, 1) - 2*leftpoints(idx, 1) - 2*rightpoints(idx+1, 1) + 2*rightpoints(idx, 1)) - 1.0/2.0*(alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))*(2*leftpoints(idx+1, 2) - 2*leftpoints(idx, 2) - 2*rightpoints(idx+1, 2) + 2*rightpoints(idx, 2)))*(atan((alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))/(alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))) - atan((-alpha(idx)*(leftpoints(idx-1, 2) - rightpoints(idx-1, 2)) + alpha(idx)*(leftpoints(idx, 2) - rightpoints(idx, 2)) - rightpoints(idx-1, 2) + rightpoints(idx, 2))/(-alpha(idx)*(leftpoints(idx-1, 1) - rightpoints(idx-1, 1)) + alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) - rightpoints(idx-1, 1) + rightpoints(idx, 1))))/((alpha(idx)*(leftpoints(idx+1, 1) - rightpoints(idx+1, 1)) - alpha(idx)*(leftpoints(idx, 1) - rightpoints(idx, 1)) + rightpoints(idx+1, 1) - rightpoints(idx, 1))^2 + (alpha(idx)*(leftpoints(idx+1, 2) - rightpoints(idx+1, 2)) - alpha(idx)*((leftpoints(idx, 2) - rightpoints(idx, 2)) + rightpoints(idx+1, 2) - rightpoints(idx, 2))^2))^(3.0/2.0);
   jr = jr+dsd;
end
J=kl*Jld+kr*Jr;

