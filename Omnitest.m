clc;
clear;
clf;

global C
cas = 2;

if cas == 1 
    C = readmatrix("Robot1.csv");
else
    M = readmatrix("Robot.csv");
    C = M(:,2);
end

t_span = linspace(0,18,100);

%%
rhoL = zeros(length(t_span),2);
rhoU = zeros(length(t_span),2);
for i=1:length(t_span)
    [rL, rU] = tube(t_span(i));
    rhoL(i,:) = rL';
    rhoU(i,:) = rU';
end

t = t_span;
%% Plots
figure(1)
dim = 1;
subplot(2,1,1)
plot(t,rhoL(:,dim),'Linestyle','-.','Color',[0.1 0.1 0.9],'Linewidth',2); hold on;
plot(t,rhoU(:,dim),'-k','Linewidth',1.5);
legend({'Low','Up'},'Fontsize',15,'Location','best')
xlabel('time','Fontsize',15)
ylabel('$x_1$ (rad)','interpreter','Latex','Fontsize',15,'Fontweight','bold')
grid on;
ax = gca;
ax.FontSize = 16;
% ylim([0 5])

dim = 2;
subplot(2,1,2)
plot(t,rhoL(:,dim),'Linestyle','-.','Color',[0.1 0.1 0.9],'Linewidth',2); hold on;
plot(t,rhoU(:,dim),'-k','Linewidth',1.5);
legend({'Low','Up'},'Fontsize',15,'Location','best')
xlabel('time','Fontsize',15)
ylabel('$x_2$ (rad)','interpreter','Latex','Fontsize',15,'Fontweight','bold')
grid on;
ax = gca;
ax.FontSize = 16;
% ylim([0 5])

%% Spatiotemporal Tubes
function [gamL, gamU] = tube(t)
    
    global C
    Cpiece = C; 

    % If the coefficients are stored as a column vector, reshape
    if size(Cpiece, 2) == 1 && mod(length(Cpiece), 4) == 0
        Cpiece = reshape(Cpiece, [], 4)';  % Reshape to 4 x N
    end

    % Piecewise breakpoints and polynomial degree
    breaks = [0 3.5 7 10.5 14 18];
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
                    gamma_Ux(k) = val;
                case 2
                    gamma_Lx(k) = val;
                case 3
                    gamma_Uy(k) = val;
                case 4
                    gamma_Ly(k) = val;
            end
        end
    end
    gamL = [gamma_Lx; gamma_Ly];
    gamU = [gamma_Ux; gamma_Uy];
end

