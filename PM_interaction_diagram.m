clear global;
clc;

% Material Properties
fck = 30; % Grade of concrete in N/mm^2
fy = 500; % Grade of steel in N/mm^2
Es = 2*10^5; % Modulus of elasticity of Steel

% Column Dimensions	
b = 500; % Width of the column in mm
D = 500; % Depth of the column in mm
dd = 50; % Effective cover depth
n = 4; % number of layers of rebars	
pt = (n-1)*pi*20*20*100/(b*D); % percentage of steel

% Rebar positions and areas
d = D - dd;
y_pos = linspace(dd, D-dd, n);
y = D/2 - y_pos; % Distance from centroidal axis
e_c = 0.0035;
e_st = -(0.002 + 0.87*fy/Es);

As = linspace(0,0,n);
for j = 1:n
    if j == 1 || j == n
        As(j) = (n)/(4*n-4)*pt*b*D/100;
    else 
        As(j) = (2)/(4*n-4)*pt*b*D/100;
    end
end
% disp(As);

xBalanced = solveNeutralAxisBalanced(e_c, e_st, d);
disp(xBalanced);
xPureMoment = solveNeutralAxisPureMoment(fy,fck,Es,y,As,D,b,d);
disp(xPureMoment);

% % % Strain conditions for different points % % %
% % Pure Compression
% e_top = 0.0035;
% e_bottom = 0.0035;
% 
% % Decompression point
% e_top = 0.0035;
% e_bottom = 0;
% 
% % Balanced point
% e_top = 0.0035;
% e_bottom = -(0.002 + 0.87*fy/Es);
% 
% % Pure Moment
% e_top = find by iteration;
% e_bottom = -(0.002 + 0.87*fy/Es);
% 
% % Pure Tension
% e_top = -(0.002 + 0.87*fy/Es);
% e_bottom = -(0.002 + 0.87*fy/Es);


% Calculate P and M for different strain conditions
points = zeros(22, 2);

i = 1;
% Pure compression
e_c = 0.0035;
fcc = fc(e_c, fck);
fsc = fs(e_c, fy);
Fc = fcc*b*D;
Fs = (fsc - fcc)*pt*b*D/100;

P = Fc + Fs;
M = 0;
points(i,1) = P/1e3;
points(i,2) = M/1e6;
i = i+1;

% Decompression point to Balanced point
xu = linspace(D,xBalanced,10);
e_c = 0.0035;
for j = 1:10
    strain = e_c/xu(j) * (y - (D/2 - xu(j)));
    % disp(strain);
    stress = arrayfun(@(e) (fs(e, fy) - fc(e,fck)) , strain);
    % disp(stress);
    fcc = fcavg(e_c,fck);
    Fc = fcc*b*xu(j);
    Fs = sum(stress .* As);
    
    P = Fc + Fs;
    M = Fc*(D/2 - xcavg(e_c, xu(j))) + sum(stress .* As .* y);
    points(i,1) = P/1e3;
    points(i,2) = M/1e6;
    i = i+1;
end

% Balanced point to Pure Moment point
xu = linspace(xBalanced,xPureMoment,10);
e_st = -(0.002 + 0.87*fy/Es);
for j = 1:10
    e_c = round(-e_st/(d-xu(j))*xu(j), 4);
    strain = e_c/xu(j) * (y - (D/2 - xu(j)));
    % disp(strain);
    stress = arrayfun(@(e) (fs(e, fy) - fc(e,fck)) , strain);
    % disp(stress);
    fcc = fcavg(e_c,fck);
    Fc = fcc*b*xu(j);
    Fs = sum(stress .* As);
    
    P = Fc + Fs;
    M = Fc*(D/2 - xcavg(e_c, xu(j))) + sum(stress .* As .* y);
    points(i,1) = P/1e3;
    points(i,2) = M/1e6;
    i = i+1;
end

% Pure Tension
P = fs(e_st, fy)*pt*b*D/100;
M = 0;
points(i,1) = P/1e3;
points(i,2) = M/1e6;
i = i+1;

% Plot P-M interaction diagram
PM_diagram(points);



% % Balanced point
% e_c = 0.0035;
% e_st = -(0.002 + 0.87*fy/Es);
% xu = xBalanced;
% strain = e_c/xu * (y - (D/2 - xu));
% % disp(strain);
% stress = arrayfun(@(e) (fs(e, fy) - fc(e,fck)) , strain);
% % disp(stress);
% fcc = fcavg(e_c,fck);
% Fc = fcc*b*xu;
% Fs = sum(stress .* As);
% 
% P = (Fc + Fs)/1e3;
% M = (Fc*(D/2 - xcavg(e_c, xu)) + sum(stress .* As .* y))/1e6;




% % Function Definitions

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

% Average Stress-strain relationship for concrete
function f = fcavg(e, fck)
    if e > 0 && e <= 0.002
        f = 0.45*fck*((e/0.002) - 1/3*(e/0.002)^2);
    elseif e > 0.002 && e <= 0.0035
        f = 0.45*fck*(1-1/3*(0.002/e));
    else 
        f = 0;
    end
end

% Depth of com of Average stress for concrete
function x = xcavg(e,X)
    if e > 0 && e <= 0.002
        x = (1/3 - ((e/0.002)/12))*X/(1-((e/0.002)/3));
    elseif e > 0.002 && e <= 0.0035
        x = (1/2 - 1/3*(0.002/e) + 1/12*((0.002/e)^2))*X/(1-((0.002/e)/3));
    else 
        x = 0;
    end
end

% Function to calculate depth of Neutral axis
function x = solveNeutralAxisBalanced(e_c, e_st, d)
    x = e_c * d/(e_c - e_st);
end

function x = solveNeutralAxisPureMoment(fy,fck,Es,y,As,D,b,d)
    % Depth of Neutral axis at Pure Moment
    e_st = -(0.002 + 0.87*fy/Es);
    e_c = 0.0035;

    x1 = round(solveNeutralAxisBalanced(e_c, e_st, d),1);
    x2 = 0;
    iter = 0;
    while(x1 ~= x2) && iter < 100
        e_c = round(-e_st/(d-x1)*x1, 4);
        strain = e_c/x1 * (y - (D/2 - x1));
        % disp(strain);
        stress = arrayfun(@(e) (fs(e, fy) - fc(e,fck)) , strain);
        % disp(stress);
        fcc = fcavg(e_c,fck);
        if fcc == 0
            error("Check, fcavg is zero");
        end
        T = -sum(stress .* As);
        x2 = T/(fcc*b);
        % x1 = round((x1 + x2)/2,1);
        x1 = x1 - 1;
        % disp(x2);
        % disp(x1);
        if abs(x2-x1) <= 2
            break
        end
        iter = iter + 1;
    end
    x = x1;
end

% Plotting PM_interaction_diagram
function PM_diagram(points)
    figure;
    plot(points(:,2), points(:,1), 'o-', 'LineWidth', 2);
    xlabel('M (kNm)');
    ylabel('P (kN)');
    title('P-M Interaction Diagram for Rectangular Column');
    grid on;
    legend('Interaction Points', 'Location', 'Best');
end

