clc;
clear;
clf;

M = readmatrix("Robot.csv");
C = M(:,2);  % 120x1 vector
breaks = [0 3 6 9 12 18];  % 6x1 vector
numSegments = length(breaks) - 1;

% Prepare to store full trajectories
Lx_all = [];
Ly_all = [];
Ux_all = [];
Uy_all = [];
t_all  = [];

for i = 1:numSegments
    % Get coefficients for this segment
    idx_start = (i-1)*24 + 1;
    coeffs = C(idx_start:idx_start+23);
    
    % Coefficients for each polynomial
    Lx_coeffs = coeffs(1:6);
    Ly_coeffs = coeffs(7:12);
    Ux_coeffs = coeffs(13:18);
    Uy_coeffs = coeffs(19:24);
    
    % Time vector for this segment
    t_local = linspace(0, breaks(i+1) - breaks(i), 100);
    t_global = t_local + breaks(i);
    
    % Evaluate polynomials
    T = [ones(size(t_local));
         t_local;
         t_local.^2;
         t_local.^3;
         t_local.^4;
         t_local.^5];

    Lx = Lx_coeffs' * T;
    Ly = Ly_coeffs' * T;
    Ux = Ux_coeffs' * T;
    Uy = Uy_coeffs' * T;

    % Store results
    Lx_all = [Lx_all, Lx];
    Ly_all = [Ly_all, Ly];
    Ux_all = [Ux_all, Ux];
    Uy_all = [Uy_all, Uy];
    t_all  = [t_all, t_global];
end

% Plotting
figure(1);
subplot(2,1,1);
plot(t_all, Lx_all, 'b', 'LineWidth', 1.5); hold on;
plot(t_all, Ux_all, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('X');
legend('Lx', 'Ux');
title('X Coordinates');

subplot(2,1,2);
plot(t_all, Ly_all, 'b', 'LineWidth', 1.5); hold on;
plot(t_all, Uy_all, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Y');
legend('Ly', 'Uy');
title('Y Coordinates');
