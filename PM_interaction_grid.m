clear global
clc

% Material Properties
fck = 30; % Grade of concrete in N/mm^2
fy = 500; % Grade of steel in N/mm^2
Es = 2*10^5; % Modulus of elasticity of Steel

% Limiting strain
e_cp = 0.0035;
e_st = 0.002+(0.87*fy/Es);

% Column Dimensions	
b = 500; % Width of the column in mm
D = 500; % Depth of the column in mm
dd = 50; % Effective cover depth
n = 4; % number of layers of rebars	
dia = 20; % bar diameter in mm

% Basic Calculations
d = D-dd;

% Steel area and position
Ast = pi*dia*dia/4;
Sx = (b-2*dd)/(n-1);
Sy = (D-2*dd)/(n-1);
st_x_pos = [-D/2+dd, -D/2+dd+Sx, -D/2+dd+2*Sx, -D/2+dd+3*Sx, -D/2+dd, -D/2+dd+3*Sx, -D/2+dd, -D/2+dd+3*Sx, -D/2+dd, -D/2+dd+Sx, -D/2+dd+2*Sx, -D/2+dd+3*Sx];
st_y_pos = [-D/2+dd, -D/2+dd, -D/2+dd, -D/2+dd, -D/2+dd+Sy, -D/2+dd+Sy, -D/2+dd+2*Sy, -D/2+dd+2*Sy, -D/2+dd+3*Sy, -D/2+dd+3*Sy, -D/2+dd+3*Sy, -D/2+dd+3*Sy];

c = 1; % dimension of one grid block

% Create the grid points
x = linspace(-b/2, b/2, b/c+1);
y = linspace(-D/2, D/2, D/c+1);
[X, Y] = meshgrid(x, y);    
% % Visualize the grid
figure;
plot(X, Y, 'bs');  % Plot grid points as blue circles
hold on;
plot(X', Y', 'r.');  % Plot transposed grid to connect the points
xlabel('X');
ylabel('Y');
title('2D Grid in a Rectangular Section');
axis equal;
grid on;
% 
% % Define the rotation angle (in degrees)
% theta = 20;
% theta_rad = deg2rad(theta);
% 
% % Rotation matrix
% R = [cos(theta_rad) -sin(theta_rad); 
%      sin(theta_rad) cos(theta_rad)];
% 
% % Apply rotation
% XY_rotated = R * [X(:)'; Y(:)'];
% 
% % Reshape rotated coordinates back into grid format
% X_rot = reshape(XY_rotated(1,:), size(X));
% Y_rot = reshape(XY_rotated(2,:), size(Y));
% 
% % Visualize the rotated grid
% figure;
% plot(X_rot, Y_rot, 'bs');  % Plot grid points as blue circles
% hold on;
% plot(X_rot', Y_rot', 'r.');  % Connect points to form the grid
% 
% axis equal;
% xlabel('X');
% ylabel('Y');
% title('Rotated Rectangular Section with 2D Grid');
% grid on;


x_pos = -b/2:c:b/2;
y_pos = -D/2:c:D/2;

y_bal = d*e_cp/(e_cp+e_st);
arr = zeros((4*D/c)+1,2); % Array to store P,M (in kN, kN-m)
count = 1;

% Varying depth of neutral axis
for yu = -D:c:3*D
    P = 0;
    M = 0;
    % Strain calculation
    if yu > y_bal
        strain = -(y_pos-(yu-D/2))*(e_cp/yu);
        strain_steel = -(st_y_pos-(yu-D/2))*(e_cp/yu);
    elseif yu <= y_bal
        strain = -(y_pos-(yu-D/2))*e_st/(d-yu);
        strain_steel = -(st_y_pos-(yu-D/2))*e_st/(d-yu);
    end

    % Stress calculation
    stress = strain;
    for i = 1:length(strain)
        stress(i) = fc(strain(i), fck);
    end
    stress_steel = strain_steel;
    for i = 1:length(strain_steel)
        stress_steel(i) = fs(strain_steel(i), fy);
    end
    
    % Concrete stress resultants
    for i = 1:(b/c)+1
        for j = 1:(D/c)+1
            P = P+stress(j)*(c*c);
            M = M+stress(j)*(c*c)*(-y_pos(j));
        end
    end

    % Steel stress resultants
    for i = 1:length(stress_steel)
        P = P+stress_steel(i)*(pi*dia*dia/4);
        M = M+stress_steel(i)*(pi*dia*dia/4)*(-st_y_pos(i));
    end
    % if count > 1 && arr(count-1,2) == 0 && M/1e6 == 0
    %     if P < 0
    %         continue;
    %     elseif P > 0
    %         break;
    %     end
    % end
    arr(count,1) = P/1e3;
    arr(count,2) = M/1e6;
    count = count + 1;
end

% Plot P-M interaction diagram
PM_diagram(arr);







% Material constitutive relations
% Stress-strain relationship for concrete
function f = fc(e,fck)
    if e > 0 && e <= 0.002
        f = 0.45*fck*(2*(e/0.002) - (e/0.002)^2);
    elseif e > 0.002 && e <= 0.0035
        f = 0.45*fck;
    else 
        f = 0;
    end
end

% Stress-strain relationship for steel
function f = fs(e,fy)
    Es = 2*10^5;
    if e >= 0 && e <= 0.87*fy/Es
        f = Es*e;
    elseif e > 0.87*fy/Es
        f = 0.87*fy;
    elseif e < 0 && e >= -0.87*fy/Es
        f = Es*e;
    elseif e < -0.87*fy/Es
        f = -0.87*fy;
    else 
        f = 0;
    end
end

% Plotting PM_interaction_diagram
function PM_diagram(arr)
    figure;
    plot(arr(:,2), arr(:,1), 'o-', 'LineWidth', 2);
    xlabel('M (kNm)');
    ylabel('P (kN)');
    title('P-M Interaction Diagram for Rectangular Column');
    grid on;
    legend('Interaction Points', 'Location', 'Best');
end