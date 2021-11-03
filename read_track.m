clear all
close all
clc
track = readtable("tracks/Monza.csv");
n = height(track);
%leftpoints = zeros(n, 2);
%rightpoints =  zeros(n, 2);
for point = 1:n-1
    %A(point, :) = get_border(track.x_X_m(point), track.y_m(point), track.x_X_m(point+1), track.y_m(point+1), track.w_tr_left_m(point), track.w_tr_left_m(point+1), track.w_tr_right_m(point), track.w_tr_right_m(point+1));
    [leftpoints(point, :), rightpoints(point, :)] = get_border(track.x_X_m(point), track.y_m(point), track.x_X_m(point+1), track.y_m(point+1), track.w_tr_left_m(point), track.w_tr_left_m(point+1), track.w_tr_right_m(point), track.w_tr_right_m(point+1));
end
hold on
plot(leftpoints(:, 1), leftpoints(:, 2))
plot(rightpoints(:, 1), rightpoints(:,2))
hold off