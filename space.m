clc;
clear;
clf;

global cas

for cas = 1:2

t_span = linspace(0,15,150);

%%
rhoL = zeros(length(t_span),3);
rhoU = zeros(length(t_span),3);
for i=1:length(t_span)
    [rL, rU] = tube(t_span(i));
    rhoL(i,:) = rL';
    rhoU(i,:) = rU';
end

X_init = 0.5*(rhoL(1,:)+rhoU(1,:))';
% X = 0.5*(rhoL+rhoU);
% t = t_span;
options = odeset('RelTol',1e-3);
tic
[t, X] = ode45(@spacecraft, t_span, X_init, options);
toc

%% Robustness
T = [2.2, 1.6, 1.6];
G = [2.8, 2.8, 2.8];
d = 0.2;
rob_arr = [];
for i = 1:length(t)
    rob = max(d - norm(X(i,:)-T), d - norm(X(i,:)-G));
    rob_arr = [rob_arr, rob];
end

%% Plots
figure(1)
dim = 1;
subplot(2,4,2 + (cas-1)*4)
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
subplot(2,4,3 + (cas-1)*4)
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
subplot(2,4,4 + (cas-1)*4)
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

subplot(2,4,1 + (cas-1)*4)
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
xlim([0,4])
ylim([0,4])
zlim([0,3])
ax = gca;
ax.FontSize = 16;
view([-15,11])

end
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
global cas
C = C_val(cas);
s = 1;
for i=1:6
    c0 = C(6*i-5);
    c1 = C(6*i-4);
    c2 = C(6*i-3);
    c3 = C(6*i-2);
    c4 = C(6*i-1);
    c5 = C(6*i);
    if i==1
        gamma_Lx = s + c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==2
        gamma_Ly = s + c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==3
        gamma_Lz = s + c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==4
        gamma_Ux = s + c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    elseif i==5
        gamma_Uy = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    else
        gamma_Uz = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
    end
end
gamL = [gamma_Lx; gamma_Ly; gamma_Lz]/5;
gamU = [gamma_Ux; gamma_Uy; gamma_Uz]/5;
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
    k = 10^3;
    u = -k*Xi*eps;
end

%% 2-Norm
function f = normTwo(a,b) 
    x = a(1)-b(1);
    y = a(2)-b(2);
    f = sqrt(x^2+y^2);
end

%%
function C = C_val(cas)
if cas == 1 
    C0 = -0.5092201124440394;
    C1 = -1.665740011109661;
    C2 = 1.542484516723121;
    C3 = -0.3243109601139466;
    C4 = 0.026244366716944464;
    C5 = -0.0007154516946045534;
    C6 = -0.5092201124440394;
    C7 = -0.9354351139789893;
    C8 = 0.8798053580061804;
    C9 = -0.15421328533867962;
    C10 = 0.011825144550615181;
    C11 = -0.00032995786104859844;
    C12 = 1.4907798875559606;
    C13 = -0.14819509743468093;
    C14 = 0.32437796115886813;
    C15 = -0.058556349964849215;
    C16 = 0.0052465859553120706;
    C17 = -0.0001699023672777245;
    C18 = 1.9953899437779803;
    C19 = -1.6883828104673564;
    C20 = 1.5955854133448955;
    C21 = -0.34080413875938786;
    C22 = 0.02799927661790214;
    C23 = -0.0007749401658234577;
    C24 = 1.9953899437779803;
    C25 = -0.9354351139789893;
    C26 = 0.8798053580061804;
    C27 = -0.15421328533867962;
    C28 = 0.011825144550615181;
    C29 = -0.00032995786104859844;
    C30 = 3.99538994377798;
    C31 = -0.14819509743468093;
    C32 = 0.32437796115886813;
    C33 = -0.058556349964849215;
    C34 = 0.0052465859553120706;
    C35 = -0.0001699023672777245;
else
    C0 = -0.9947193621508238;
    C1 = 0.4277863285546202;
    C2 = 1.2140297298487956;
    C3 = -0.27679682188266586;
    C4 = 0.02200792560735738;
    C5 = -0.0005888928692527835;
    C6 = -0.5105612756983524;
    C7 = -1.4667281217292212;
    C8 = 1.0843122115171406;
    C9 = -0.18152804910001696;
    C10 = 0.013280466284529367;
    C11 = -0.0003552471176615102;
    C12 = 1.4894387243016476;
    C13 = -1.27618875247245;
    C14 = 0.655150320834794;
    C15 = -0.08099965733727106;
    C16 = 0.004422992408222249;
    C17 = -8.864906559813607e-05;
    C18 = 1.5105612756983524;
    C19 = 0.4277863285546202;
    C20 = 1.2140297298487956;
    C21 = -0.27679682188266586;
    C22 = 0.02200792560735738;
    C23 = -0.0005888928692527835;
    C24 = 1.994719362150824;
    C25 = -1.4667281217292212;
    C26 = 1.0843122115171406;
    C27 = -0.18152804910001696;
    C28 = 0.013280466284529367;
    C29 = -0.0003552471176615102;
    C30 = 3.994719362150824;
    C31 = -1.27618875247245;
    C32 = 0.655150320834794;
    C33 = -0.08099965733727106;
    C34 = 0.004422992408222249;
    C35 = -8.864906559813607e-05;
end
C = [C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, ...
        C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26, C27, C28, C29, C30, C31, C32, C33, C34, C35];
end