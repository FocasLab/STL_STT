clc;
clear;
clf;

global C
cas = 2;

if cas == 1 
    C = readmatrix("Robot1.csv");
else
    C = readmatrix("Robot.csv");
end

t_span = linspace(0.4,18,1000);

%%
rhoL = zeros(length(t_span),2);
rhoU = zeros(length(t_span),2);
for i=1:length(t_span)
    [rL, rU] = tube(t_span(i));
    rhoL(i,:) = rL';
    rhoU(i,:) = rU';
end

X_init = 0.5*(rhoL(1,:)+rhoU(1,:))';
options = odeset('RelTol',1e-5);

tic
[t, X] = ode45(@spacecraft, t_span, X_init, options);
toc

% X = 0.5*(rhoL+rhoU);
% t = t_span;

%% Plots
figure(1)
dim = 1;
subplot(2,2,2)
plot(t,rhoL(:,dim),t,rhoU(:,dim),'Linestyle','-.','Color',[0.1 0.1 0.9],'Linewidth',2); hold on;
plot(t,X(:,dim),'-k','Linewidth',1.5);
legend({'STTs','','System Trajectory'},'Fontsize',15,'Location','best')
xlabel('time','Fontsize',15)
ylabel('$x_1$ (rad)','interpreter','Latex','Fontsize',15,'Fontweight','bold')
grid on;
ax = gca;
ax.FontSize = 16;
% ylim([0 5])

dim = 2;
subplot(2,2,4)
plot(t,rhoL(:,dim),t,rhoU(:,dim),'Linestyle','-.','Color',[0.1 0.1 0.9],'Linewidth',2); hold on;
plot(t,X(:,dim),'-k','Linewidth',1.5);
legend({'STTs','','System Trajectory'},'Fontsize',15,'Location','best')
xlabel('time','Fontsize',15)
ylabel('$x_2$ (rad)','interpreter','Latex','Fontsize',15,'Fontweight','bold')
grid on;
ax = gca;
ax.FontSize = 16;
% ylim([0 5])

subplot(2,2,[1,3])
% for i=1:1000
    % t(i);
    i = 1000;
hold on;
rectangle('Position', [4 4 1 1], 'FaceColor', [0,0,1,0.5],'EdgeColor','none', FaceAlpha=0.5) % Start
rectangle('Position', [8 14 1 1], 'FaceColor', [0,1,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Target1
rectangle('Position', [14 8 1 1], 'FaceColor', [0,1,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Target2
rectangle('Position', [5 8 2 2], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.6) % Obstacle1
rectangle('Position', [15 12 2 2], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.6) % Obstacle2
rectangle('Position', [18 18 1 1], 'FaceColor', [0,1,0,0.5],'EdgeColor','none', FaceAlpha=0.8) % Goal
plot(X(:,1),X(:,2),'k-','Linewidth',1.5);
xlabel('$x_1$ (m)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
ylabel('$x_2$ (m)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
grid on;
box on;
xlim([0,20])
ylim([0,20])
ax = gca;
ax.FontSize = 16;
axis square;

%% System Dynamics
function dXdt = spacecraft(t, X)
    u = real(control(t,X));
    dx1 = u(1);
    dx2 = u(2);
    dXdt = [dx1; dx2];
end

%% Control Law
function u = control(t,X)   
    % Funnel
    [rhoL, rhoU] = tube(t); % interpolate value at t
    rho_a = rhoL + rhoU;
    rho_m = -rhoL + rhoU;

    % Normalized Error
    e = (X-0.5*rho_a)./(0.5*rho_m);
    
    % Transformed Error
    eps = log((1+e)./(1-e));
    
    % Xi
    Xi_col = 4./rho_m./(1-e.^2);
    Xi = diag(Xi_col);
    
    % Control Law
    k = 1;
    u = -k*Xi*eps;
end

%% Spatiotemporal Tubes
function [gamL, gamU] = tube(t)
    
    global C
    Cpiece = C; 

    % If the coefficients are stored as a column vector, reshape
    if size(Cpiece, 2) == 1 && mod(length(Cpiece), 4) == 0
        Cpiece = reshape(Cpiece, [], 4)';  % Reshape to 4 x N
    end

    % Piecewise breakpoints and polynomial degree
    breaks = [0 4 7 10 13 18];
    num_segments = length(breaks) - 1;
    degree = 5;
    coefs_per_seg = degree + 1;

    % Ensure t is a column vector
    t = t(:);

    % Initialize tube boundaries
    gamma_Lx = zeros(size(t));
    gamma_Ly = zeros(size(t));
    gamma_Ux = zeros(size(t));
    gamma_Uy = zeros(size(t));

    for k = 1:length(t)
        tk = t(k);

        % Clamp t to valid range
        if tk < breaks(1)
            tk = breaks(1);
        elseif tk >= breaks(end)
            tk = breaks(end) - 1e-8;
        end

        % Find segment
        seg_idx = find(breaks <= tk, 1, 'last');
        if seg_idx >= num_segments + 1
            seg_idx = num_segments;
        end

        t_local = tk - breaks(seg_idx);

        % Evaluate each curve
        for curve = 1:4
            coefs = Cpiece(curve, (seg_idx-1)*coefs_per_seg + (1:coefs_per_seg));
            val = polyval(flip(coefs), t_local);
            switch curve
                case 1
                    gamma_Lx(k) = val;
                case 2
                    gamma_Ux(k) = val;
                case 3
                    gamma_Ly(k) = val;
                case 4
                    gamma_Uy(k) = val;
            end
        end
    end
    gamL = [gamma_Lx; gamma_Ly];
    gamU = [gamma_Ux; gamma_Uy];
end

