  function [J, grad] = costfunction_laptime(n,rightpoints,leftpoints,alpha,axmax, aymax)

% MYCOSTFUNCTION computes the value of the cost function. 
%  INPUTS:  n             = number of segments dividing the track
%           rightpoints   = coordinates of the inner points of the track
%           leftpoints    = coordinates of the outer points of the track
%           alpha         = optimization variable, needed to express where
%                          the car has to stay inside the track
%           gamma    = weight if out of constraints
%
%  OUTPUTS: J             = value of the cost function
%               gradient     = value of the cost function gradient
[x,y] = mycoordinates(n,rightpoints,leftpoints,alpha);
J = laptime(x, y, n, axmax, aymax);
grad = 0;



