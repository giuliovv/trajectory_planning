clear all
close all 
clc
% Enter coordinates of points A and B and  
% the 1/2 length of the perpendicular segment
track = readtable("tracks/Monza.csv");
A = [track.x_X_m(1) track.y_m(1)]; %[x,y]
B = [track.x_X_m(2),track.y_m(2)]; %[x,y]
Clen = (track.w_tr_left_m(1)+track.w_tr_left_m(2) + track.w_tr_right_m(1)+track.w_tr_right_m(2))/2;  %length of line CD (1/2 of the full perpendicular line)
%% Do the math
% Get slope and y int of line AB
slope = (B(2)-A(2)) / (B(1)-A(1)); 
yint = B(2) - slope*B(1); 
% Choose a point C along line AB at half distance
C(1) = range([B(1),A(1)])/2+min([A(1),B(1)]); 
C(2) = slope * C*(1) + yint; 
% Get slope and y int of line perpendicular to AB at point C
perpSlope = -1/slope; 
perpYint = C(2) - perpSlope*C(1); 
% Find the end points of the perpendicular line with length Clen*2
x = C(1) + (Clen*sqrt(1/(1+perpSlope^2)))*[-1,1]; 
y = C(2) + (perpSlope*Clen*sqrt(1/(1+perpSlope^2)))*[-1,1]; 
%% Plot results
figure()
% Draw line AB
p1 = plot([A(1),B(1)], [A(2),B(2)], 'k-o','LineWidth',2,'DisplayName','AB'); 
hold on
axis equal 
grid on
axlim = xlim(gca)+[-1,1]*range(xlim(gca))*1.1; 
xlim(axlim)
ylim(axlim)
% Draw a references line parallel to AB
rl = refline(slope,yint); 
set(rl, 'LineStyle', '--', 'Color', [.5 .5 .5])
% Draw a references line perpendicular to AB at point C
rl = refline(perpSlope,perpYint); 
set(rl, 'LineStyle', '--', 'Color', [.5 .5 .5])
xlim(axlim) %refline() changes limits so we set them again :(
ylim(axlim)
% Add point C
p2 = plot(C(1),C(2), 'r*','DisplayName','C');
% Add the perp segment
p3 = plot(x,y,'b-o','LineWidth',2,'DisplayName','perp AB');
legend([p1,p2,p3])