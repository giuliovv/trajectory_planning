clear all
close all
clc
track = readtable("tracks/Monza.csv");
n = height(track);
leftpoints = zeros(n);
rightpoints =  zeros(n);
for point = 1:n-1
    [leftpoints(point), rightpoints(point)] = get_border(track.x_X_m(point), track.y_m(point), track.x_X_m(point+1), track.y_m(point+1), track.w_tr_left_m(point), track.w_tr_left_m(point+1), track.w_tr_right_m(point), track.w_tr_right_m(point+1));
end
hold on
plot(leftpoints)
plot(rightpoints)
hold off