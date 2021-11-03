clear all
close all
clc
track = readtable("tracks/Sochi.csv");
n = height(track);
leftpoints = zeros(n-1, 2);
rightpoints =  zeros(n-1, 2);
oldleft = [0,0];
oldright = [0,0];
for idx = 1:n-1
    [leftpoints(idx, :), rightpoints(idx, :)] = get_border(track.x_X_m(idx), track.y_m(idx), track.x_X_m(idx+1), track.y_m(idx+1), track.w_tr_left_m(idx), track.w_tr_left_m(idx+1), track.w_tr_right_m(idx), track.w_tr_right_m(idx+1), oldleft, oldright);
    oldleft(:) = leftpoints(idx, :);
    oldright(:) = rightpoints(idx, :);
end

hold on
plot(leftpoints(:, 1), leftpoints(:, 2))
plot(rightpoints(:, 1), rightpoints(:,2))
hold off