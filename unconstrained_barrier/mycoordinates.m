function [x,y] = mycoordinates(n,rightpoints,leftpoints,alpha)

% MYCOORDINATES computes the coordinates x and y of the points through which the car has to pass at each step. 
%  INPUTS:  n             = number of segments dividing the track
%           rightpoints   = coordinates of the inner points of the track
%           leftpoints    = coordinates of the outer points of the track
%           alpha         = optimization variable, needed to express where
%                          the car has to stay inside the track
%
%  OUTPUTS: x             = values of x at each step
%           y             = values of y at each step

for idx = 1:n-1
    x(idx)=rightpoints(idx,1)+alpha(idx)*(leftpoints(idx,1)-rightpoints(idx,1));
    y(idx)=rightpoints(idx,2)+alpha(idx)*(leftpoints(idx,2)-rightpoints(idx,2));
end   
