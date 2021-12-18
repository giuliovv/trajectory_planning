clear all
close all

%% True curve:
center=[10;5];
radius=7;
theta_var=linspace(0,2*pi,1e3);
points=[center(1,1)+cos(theta_var)*radius;center(2,1)+sin(theta_var)*radius];

%% Test with random points
theta=sort(rand(1,3))*2*pi;
data_points=[center(1,1)+cos(theta)*radius;center(2,1)+sin(theta)*radius];

figure,plot(points(1,:),points(2,:),data_points(1,:),data_points(2,:),'*',center(1,1),center(2,1),'*r')

A=[data_points' ones(3,1)]+eye(3)*1e-15;
b=sum(data_points'.^2,2);

param=A\b;
center_hat=[param(1,1)/2;param(2,1)/2]
radius_hat=sqrt(sum(center_hat.^2)+param(3,1))
curvature_hat=1/radius_hat

%% Equation based on arctangent
curvature_2=2*abs(atan2(data_points(2,end)-data_points(2,2),data_points(1,end)-data_points(1,2))...
-atan2(data_points(2,2)-data_points(2,1),data_points(1,2)-data_points(1,1)))...
/(norm(data_points(:,end)-data_points(:,1)))


%% Test with close points

theta=pi/4+sort(rand(1,3)*0.01);    % points are maximum 0.01 radiants apart, with a minimum of pi/4
data_points=[center(1,1)+cos(theta)*radius;center(2,1)+sin(theta)*radius];

figure,plot(points(1,:),points(2,:),data_points(1,:),data_points(2,:),'-*',center(1,1),center(2,1),'*r')

A=[data_points' ones(3,1)]+eye(3)*1e-15;
b=sum(data_points'.^2,2);

param=A\b;
center_hat=[param(1,1)/2;param(2,1)/2]
radius_hat=sqrt(sum(center_hat.^2)+param(3,1))
curvature_hat=1/radius_hat

%% Equation based on arctangent
curvature_2=2*abs(atan2(data_points(2,end)-data_points(2,2),data_points(1,end)-data_points(1,2))...
-atan2(data_points(2,2)-data_points(2,1),data_points(1,2)-data_points(1,1)))...
/(norm(data_points(:,end)-data_points(:,1)))


%% Test with points aligned on a line

theta=sort(rand(1,2))*2*pi;    % 2 points are maximum 1 radiants apart, with a minimum of pi/4
data_points=[center(1,1)+cos(theta)*radius;center(2,1)+sin(theta)*radius];
data_points=[data_points(:,1) 0.5*(data_points(:,1)+data_points(:,2)) data_points(:,2)]; % Third point is in between the other two (on a line)

figure,plot(data_points(1,:),data_points(2,:),'-*',center(1,1),center(2,1),'*r')

A=[data_points' ones(3,1)]+eye(3)*1e-15;
b=sum(data_points'.^2,2);

param=A\b;
center_hat=[param(1,1)/2;param(2,1)/2]
radius_hat=sqrt(sum(center_hat.^2)+param(3,1))
curvature_hat=1/radius_hat

%% Equation based on arctangent
curvature_2=2*abs(atan2(data_points(2,end)-data_points(2,2),data_points(1,end)-data_points(1,2))...
-atan2(data_points(2,2)-data_points(2,1),data_points(1,2)-data_points(1,1)))...
/(norm(data_points(:,end)-data_points(:,1)))

