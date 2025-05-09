clc;
clear;
clf;

global C
cas = 2;

if cas == 1 
    M = readmatrix("Spacecraft1.csv");
    C = M(:,2);
elseif cas ==2
    M = readmatrix("Spacecraft2.csv");
    C = M(:,2);
else
    M = readmatrix("Spacecraft.csv");
    C = M(:,2);
end

t_span = linspace(0,14,150);

%%
rhoL = zeros(length(t_span),3);
rhoU = zeros(length(t_span),3);
for i=1:length(t_span)
    [rL, rU] = tube(t_span(i));
    rhoL(i,:) = rL';
    rhoU(i,:) = rU';
end

X_init = 0.5*(rhoL(1,:)+rhoU(1,:));

options = odeset('RelTol',1e-5);
tic
[t, X] = ode45(@spacecraft, t_span, X_init, options);
toc

%% Robustness
T = [2.2, 1.6, 1.6];
G = [2.8, 2.8, 2.8];
d = 0.2;

%% Plots
figure(1)
dim = 1;
subplot(1,4,2)
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
subplot(1,4,3)
plot(t,rhoL(:,dim),t,rhoU(:,dim),'Linestyle','-.','Color',[0.1 0.1 0.9],'Linewidth',2); hold on;
plot(t,X(:,dim),'-k','Linewidth',1.5);
legend({'STTs','','System Trajectory'},'Fontsize',15,'Location','best')
xlabel('time','Fontsize',15)
ylabel('$x_2$ (rad)','interpreter','Latex','Fontsize',15,'Fontweight','bold')
grid on;
ax = gca;
ax.FontSize = 16;
% ylim([0 5])

dim = 3;
subplot(1,4,4)
% plot(rob_arr)
plot(t,rhoL(:,dim),t,rhoU(:,dim),'Linestyle','-.','Color',[0.1 0.1 0.9],'Linewidth',2); hold on;
plot(t,X(:,dim),'-k','Linewidth',1.5);
legend({'STTs','','System Trajectory'},'Fontsize',15,'Location','best')
xlabel('time','Fontsize',15)
ylabel('$x_3$ (rad)','interpreter','Latex','Fontsize',15,'Fontweight','bold')
grid on;
ax = gca;
ax.FontSize = 16;
% ylim([0 5])

subplot(1,4,1)
s = 2;
hold on;
plotcube([0.6,0.6,0.6],[0,0,0],0.4,[0,0,1]); % start
plotcube([0.6,0.6,0.6],[2,1.4,1.4],0.3,[0,1,0]); % T1
plotcube([0.6,0.6,0.6],[0.8,1.4,1.4],0.3,[0,1,0]); % T2
plotcube([0.6,0.6,0.6],[2.6,2.6,2.6],0.5,[0,1,0]); % goal
plotcube([0.6,1,3],[1.4,1.4,0],0.4,[1,0,0]); % obstacle

plot3(X(:,1),X(:,2),X(:,3),'k-','Linewidth',1.5);
xlabel('$x_1$ (rad)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
ylabel('$x_2$ (rad)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
zlabel('$x_2$ (rad)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
grid on;
box on;
% xlim([0,4])
% ylim([0,4])
% zlim([0,3])
ax = gca;
ax.FontSize = 16;
view([-15,11])

%% System Dynamics
function dXdt = spacecraft(t, X)
    J1 = 200;
    J2 = 200;
    J3 = 100;
    u = real(control(t,X));
    dx1 = (J2-J3)/J1*X(2)*X(3) + 1/J1*u(1);
    dx2 = (J3-J1)/J2*X(1)*X(3) + 1/J2*u(2);
    dx3 = (J1-J2)/J3*X(2)*X(1) + 1/J3*u(3);
    dXdt = [dx1; dx2; dx3];
end

%% Spatiotemporal Tubes
function [gamL, gamU] = tube(t)
global C
% C = C_val(cas);
for i=1:6
    c0 = C(6*i-5);
    c1 = C(6*i-4);
    c2 = C(6*i-3);
    c3 = C(6*i-2);
    c4 = C(6*i-1);
    c5 = C(6*i);
    if i==1
        gamma_Lx = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==2
        gamma_Ly = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==3
        gamma_Lz = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==4
        gamma_Ux = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==5
        gamma_Uy = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    else
        gamma_Uz = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    end
end
gamL = [gamma_Lx; gamma_Ly; gamma_Lz];
gamU = [gamma_Ux; gamma_Uy; gamma_Uz];
end

%% Control Law
function u = control(t,X)   
    % STT
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
    k = 100;
    u = -k*Xi*eps;
end

%% 2-Norm
function f = normTwo(a,b) 
    x = a(1)-b(1);
    y = a(2)-b(2);
    f = sqrt(x^2+y^2);
end
