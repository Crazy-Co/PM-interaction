clear global
clc

% Input -------------------------------------------------------
% Material Properties
fck = 30; % Grade of concrete in N/mm^2
fy = 500; % Grade of steel in N/mm^2
Es = 2*10^5; % Modulus of elasticity of Steel

% Column Dimensions	
b = 400; % Width of the column in mm
D = 400; % Depth of the column in mm
dd = 50; % Effective cover depth
n = 4; % number of layers of rebars	
dia = 20; % bar diameter in mm

% Limiting strain
e_cp = 0.0035;
e_st = round(0.002+(0.87*fy/Es),4);

% Basic Calculations -----------------------------------------------------
% d = D-dd;

% theta = 30; % Define the rotation angle (in degrees)
step = 0.0001;
r = 5; % rotate neutral axis by 'r' degree
k = 1; % iterator
arr = zeros((round(e_cp+e_st,4)*2/step+1)*(90/r+1),3); % Array to store P,M (in kN, kN-m)

for theta = 0:r:90
    thetaRad = deg2rad(theta);
    
    % Rotation matrix
    R = [cos(thetaRad) sin(thetaRad); 
         -sin(thetaRad) cos(thetaRad)];
    
    % Steel area -------------------------------------------------------
    Ast = pi*dia*dia/4;
    As = zeros(n, n);  % Assuming As is initialized as a zero matrix
    
    % Assign value to the boundaries directly
    As(1, :) = Ast;     % Top row
    As(n, :) = Ast;     % Bottom row
    As(:, 1) = Ast;     % Left column
    As(:, n) = Ast;     % Right column
    
    % Steel position -------------------------------------------------------
    steel_x = linspace(-b/2+dd, b/2-dd,n);
    steel_y = linspace(-D/2+dd, D/2-dd,n);
    [Xs, Ys] = meshgrid(steel_x, steel_y);
    
    % Apply rotation
    XsYs_rotated = R * [Xs(:)'; Ys(:)'];
    % Reshape rotated coordinates back into grid format
    Xs = reshape(XsYs_rotated(1,:), size(Xs));
    Ys = reshape(XsYs_rotated(2,:), size(Ys));
    
    
    % Concrete meshing -------------------------------------------------------
    c = 5; % dimension of one grid block
    
    % Create the grid points
    x = linspace(-b/2, b/2, b/c+1);
    y = linspace(-D/2, D/2, D/c+1);
    [Xc, Yc] = meshgrid(x, y);    
    % % Visualize the grid
    % figure;
    % plot(X, Y, 'bs');  % Plot grid points as blue circles
    % hold on;
    % plot(X', Y', 'r.');  % Plot transposed grid to connect the points
    % xlabel('X');
    % ylabel('Y');
    % title('2D Grid in a Rectangular Section');
    % axis equal;
    % grid on;
    
    % Apply rotation
    XcYc_rotated = R * [Xc(:)'; Yc(:)'];
    
    % Reshape rotated coordinates back into grid format
    Xc = reshape(XcYc_rotated(1,:), size(Xc));
    Yc = reshape(XcYc_rotated(2,:), size(Yc));
    
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
    
    
    % P-M interaction calculation for given theta
    e_top = e_cp;
    
    for e_bottom = e_cp:-step:-e_st
        % P, M effect due to concrete --------------------------------------------
        % row-wise traversal
        e = zeros(size(Yc));
        f = zeros(size(Yc));
        P = 0;
        M = 0;

        for i=1:D/c+1 % for rows traversal 
            for j=1:b/c+1 % for column traversal
                e(i,j) = e_bottom + ((e_top-e_bottom)*(Yc(i,j)-Ys(1,n))/(Yc(size(Yc,1),1)-Ys(1,n)));
                f(i,j) = fc(e(i,j),fck);
                P = P + f(i,j)*(c*c);
                M = M + f(i,j)*(c*c)*Yc(i,j);
            end
        end
        
        % P, M effect due to steel -----------------------------------------------
        est = zeros(n,n);
        fst = zeros(n,n);
        for i = 1:n
            for j = 1:n
                if As(i,j) ~= 0
                    est(i,j) = e_bottom + ((e_top-e_bottom)*(Ys(i,j)-Ys(1,n))/(Yc(size(Yc,1),1)-Ys(1,n)));
                    fst(i,j) = fs(est(i,j),fy);
                    P = P + fst(i,j)*As(i,j);
                    M = M + fst(i,j)*As(i,j)*Ys(i,j);
                end
            end
        end
        fprintf("(%d, %d)\n", P/1e3, M/1e6);
        arr(k,1) = P/1e3;
        arr(k,2) = M/1e6;
        arr(k,3) = theta;
        k = k+1;
    end
    
    for e_top = e_cp:-step:-e_st
        % P, M effect due to concrete --------------------------------------------
        % row-wise traversal
        e = zeros(size(Yc));
        f = zeros(size(Yc));
        P = 0;
        M = 0;
        for i=1:D/c+1 % for rows traversal 
            for j=1:b/c+1 % for column traversal
                e(i,j) = e_bottom + ((e_top-e_bottom)*(Yc(i,j)-Ys(1,n))/(Yc(size(Yc,1),1)-Ys(1,n)));
                f(i,j) = fc(e(i,j),fck);
                P = P + f(i,j)*(c*c);
                M = M + f(i,j)*(c*c)*Yc(i,j);
            end
        end
        
        % P, M effect due to steel -----------------------------------------------
        est = zeros(n,n);
        fst = zeros(n,n);
        for i = 1:n
            for j = 1:n
                if As(i,j) ~= 0
                    est(i,j) = e_bottom + ((e_top-e_bottom)*(Ys(i,j)-Ys(1,n))/(Yc(size(Yc,1),1)-Ys(1,n)));
                    fst(i,j) = fs(est(i,j),fy);
                    P = P + fst(i,j)*As(i,j);
                    M = M + fst(i,j)*As(i,j)*Ys(i,j);
                end
            end
        end
        fprintf("(%d, %d)\n", P/1e3, M/1e6);
        arr(k,1) = P/1e3;
        arr(k,2) = M/1e6;
        arr(k,3) = theta;
        k = k+1;
    end
end    
    
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
end