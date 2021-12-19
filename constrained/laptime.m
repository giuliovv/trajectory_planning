function t = laptime(x, y, n, axmax, aymax)

rho = zeros(n,1);
ds = zeros(n,1);

v = zeros(n,1);
ay = zeros(n,1);
ax = zeros(n,1);

ds(1)=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);

% Initial guess
for idx = 2:n-2
    ds(idx)=sqrt((x(idx+1)-x(idx))^2+(y(idx+1)-y(idx))^2);
    dtheta=atan2((y(idx+1)-y(idx)),(x(idx+1)-x(idx)))-atan2((y(idx)-y(idx-1)),(x(idx)-x(idx-1)));
    rho(idx) = dtheta/ds(idx);
    v(idx) = sqrt(aymax/abs(rho(idx)));
end

v2 = v;

% Forward pass
for idx = 2:n-2
    ay(idx-1) = v2(idx-1)^2*rho(idx-1);
    ax(idx-1) = real(axmax*sqrt(1-(ay(idx-1)/aymax)^2));
    vguess = real(sqrt(v2(idx-1)^2 + 2*ax(idx-1)*ds(idx-1)));
    v2(idx) = min(v(idx), vguess);
end

v3 = v2;

% Backward pass
for idx = flip(2:n-2)
    ay(idx+1) = v3(idx+1)^2*rho(idx+1);
    ax(idx+1) = real(axmax*sqrt(1-(ay(idx+1)/aymax)^2));
    vguess = real(sqrt(v3(idx+1)^2 + 2*ax(idx+1)*ds(idx)));
    v3(idx) = min(v2(idx), vguess);
end

% Time calc
t = 0;
for idx = 1:n-1
    if ax(idx) == 0
        continue;
    end
    t = t + (-v3(idx)+sqrt(v3(idx)^2+2*ax(idx)*ds(idx)))/ax(idx);
end