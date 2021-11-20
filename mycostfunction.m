  function [J, gradient]= mycostfunction(n,rightpoints,leftpoints,alpha)

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
kr=0;


[x,y] = mycoordinates(n,rightpoints,leftpoints,alpha);

%cost function

F = zeros(n);

Jl=0;
Jld=zeros(n);
for idx = 1:n-2
    ds=(x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2;
    Jl=Jl+ds;
    
    F(idx, idx) = kl*ds;
    
    % Jacobian
    Jld(idx, idx) = (-2*leftpoints(idx, 1) + 2*rightpoints(idx, 1))*((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1)) + (-2*leftpoints(idx, 2) + 2*rightpoints(idx, 2))*((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2));
    Jld(idx, idx+1) = (2*leftpoints(idx+1, 1) - 2*rightpoints(idx+1, 1))*((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1)) + (2*leftpoints(idx+1, 2) - 2*rightpoints(idx+1, 2))*((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2));
    
end

Jr=0;
Jrd=zeros(n);
for idx = 2:n-2
    ds=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    dtheta=atan((y(idx+1)-y(idx))/(x(idx+1)-x(idx)))-atan((y(idx)-y(idx-1))/(x(idx)-x(idx-1)));
    rho=dtheta/ds;
    Jr=Jr+rho^2;
    
    F(idx+1, idx+1) = F(idx+1, idx+1)+kr*rho^2;
    
     % Jacobian
    Jrd(idx+1, idx+1) =(-((-leftpoints(idx, 1) + rightpoints(idx, 1))*(-(leftpoints(idx-1, 2) - rightpoints(idx-1, 2))*alpha(idx-1) + (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) - rightpoints(idx-1, 2) + rightpoints(idx, 2))/power(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1), 2) + (leftpoints(idx, 2) - rightpoints(idx, 2))/(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1)))/(1 + power(-(leftpoints(idx-1, 2) - rightpoints(idx-1, 2))*alpha(idx-1) + (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) - rightpoints(idx-1, 2) + rightpoints(idx, 2), 2)/power(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1), 2)) + ((leftpoints(idx, 1) - rightpoints(idx, 1))*((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2))/power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + (-leftpoints(idx, 2) + rightpoints(idx, 2))/((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1)))/(1 + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2)/power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2)))/sqrt(power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2)) + (-1.0/2.0*(-2*leftpoints(idx, 1) + 2*rightpoints(idx, 1))*((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1)) - 1.0/2.0*(-2*leftpoints(idx, 2) + 2*rightpoints(idx, 2))*((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2))) ...
        *(atan(((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2))/((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1))) - atan((-(leftpoints(idx-1, 2) - rightpoints(idx-1, 2))*alpha(idx-1) + (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) - rightpoints(idx-1, 2) + rightpoints(idx, 2))/(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1))))/power(power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2), 3.0/2.0);
    Jrd(idx+1, idx+2) = (-1.0/2.0*(2*leftpoints(idx+1, 1) - 2*rightpoints(idx+1, 1))*((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1)) - 1.0/2.0*(2*leftpoints(idx+1, 2) - 2*rightpoints(idx+1, 2))*((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2)))*(atan(((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2))/((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1))) - atan((-(leftpoints(idx-1, 2) - rightpoints(idx-1, 2))*alpha(idx-1) + (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) - rightpoints(idx-1, 2) + rightpoints(idx, 2))/(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1))))/power(power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2), 3.0/2.0) + ((-leftpoints(idx+1, 1) + rightpoints(idx+1, 1))*((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2))/power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + (leftpoints(idx+1, 2) - rightpoints(idx+1, 2))/((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1)))/((1 + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2)/power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2))*sqrt(power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2)));
    Jrd(idx+1, idx) = -((leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*(-(leftpoints(idx-1, 2) - rightpoints(idx-1, 2))*alpha(idx-1) + (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) - rightpoints(idx-1, 2) + rightpoints(idx, 2))/power(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1), 2) + (-leftpoints(idx-1, 2) + rightpoints(idx-1, 2))/(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1)))/((1 + power(-(leftpoints(idx-1, 2) - rightpoints(idx-1, 2))*alpha(idx-1) + (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) - rightpoints(idx-1, 2) + rightpoints(idx, 2), 2)/power(-(leftpoints(idx-1, 1) - rightpoints(idx-1, 1))*alpha(idx-1) + (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) - rightpoints(idx-1, 1) + rightpoints(idx, 1), 2))*sqrt(power((leftpoints(idx+1, 1) - rightpoints(idx+1, 1))*alpha(idx+1) - (leftpoints(idx, 1) - rightpoints(idx, 1))*alpha(idx) + rightpoints(idx+1, 1) - rightpoints(idx, 1), 2) + power((leftpoints(idx+1, 2) - rightpoints(idx+1, 2))*alpha(idx+1) - (leftpoints(idx, 2) - rightpoints(idx, 2))*alpha(idx) + rightpoints(idx+1, 2) - rightpoints(idx, 2), 2)));
    
end
%J=abs(kl*Jl+kr*Jr);

J=F'*F;
gradF=(kl*Jld+kr*Jrd)';
gradient = gradF*F;

%gradient = kl*Jld+kr*Jrd;

