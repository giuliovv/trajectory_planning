clear all
close all
clc

[rightpoints, leftpoints] = read_track("../tracks/Monza.csv");

hold on
plot(leftpoints(:, 1), leftpoints(:, 2), 'k')
plot(rightpoints(:, 1), rightpoints(:,2), 'k')
hold off