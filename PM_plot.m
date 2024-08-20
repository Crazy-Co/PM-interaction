% Extract P, M, theta from data
P = arr(:, 1);
M = arr(:, 2);
theta_deg = arr(:, 3);

% Convert theta from degrees to radians
theta_rad = deg2rad(theta_deg);

% Calculate Mx and My
Mx = M .* cos(theta_rad);
My = M .* sin(theta_rad);

% Plot P on z-axis, Mx on x-axis, and My on y-axis
figure;
scatter3(Mx, My, P, 5, 'filled');
xlabel('Mx');
ylabel('My');
zlabel('P');
title('3D Plot of P vs Mx and My');
grid on;