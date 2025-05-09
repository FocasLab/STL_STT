M = readmatrix("Robot.csv");
C = M(:,2);  % 120x1 vector
breaks = [0 3.5 7 10.5 14 18];
numSegments = length(breaks) - 1;

% Store full trajectories
Lx_all = [];
Ly_all = [];
Ux_all = [];
Uy_all = [];
t_all  = [];

for i = 1:numSegments
    idx_start = (i-1)*24 + 1;
    coeffs = C(idx_start:idx_start+23);

    Lx_coeffs = coeffs(1:6);
    Ly_coeffs = coeffs(7:12);
    Ux_coeffs = coeffs(13:18);
    Uy_coeffs = coeffs(19:24);

    t_local = linspace(0, breaks(i+1) - breaks(i), 100);
    t_global = t_local + breaks(i);

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

    Lx_all = [Lx_all, Lx];
    Ly_all = [Ly_all, Ly];
    Ux_all = [Ux_all, Ux];
    Uy_all = [Uy_all, Uy];
    t_all  = [t_all, t_global];
end

% Fit degree-16 polynomials
deg = 16;
pLx = polyfit(t_all, Lx_all, deg);
pLy = polyfit(t_all, Ly_all, deg);
pUx = polyfit(t_all, Ux_all, deg);
pUy = polyfit(t_all, Uy_all, deg);

% Evaluate fitted polynomials
t_fine = linspace(breaks(1), breaks(end), 1000);
Lx_fit = polyval(pLx, t_fine);
Ly_fit = polyval(pLy, t_fine);
Ux_fit = polyval(pUx, t_fine);
Uy_fit = polyval(pUy, t_fine);

% Plot fitted polynomials
figure;
subplot(2,1,1);
plot(t_fine, Lx_fit, 'b', 'LineWidth', 1.5); hold on;
plot(t_fine, Ux_fit, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('X');
legend('Lx (fit)', 'Ux (fit)');
title('Degree-16 Polynomial Fit - X Coordinates');

subplot(2,1,2);
plot(t_fine, Ly_fit, 'b', 'LineWidth', 1.5); hold on;
plot(t_fine, Uy_fit, 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Y');
legend('Ly (fit)', 'Uy (fit)');
title('Degree-16 Polynomial Fit - Y Coordinates');

% Optional: 2D envelope plot
figure;
plot(Lx_fit, Ly_fit, 'b', 'LineWidth', 1.5); hold on;
plot(Ux_fit, Uy_fit, 'r', 'LineWidth', 1.5);
xlabel('X'); ylabel('Y');
legend('Lower Bound', 'Upper Bound');
title('Envelope in XY Plane (Degree-16 Fit)');
axis equal;
