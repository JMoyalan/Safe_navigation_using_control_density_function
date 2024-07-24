%% Obstacle Avoidance control using quadratic program
clc; clear; %close all; 
mkdir('functions');
addpath('./functions');

%% Initialization
nav_p.x0 = [0;1.5];%[-1;-1.5]; % Start position
%[0;1.5]
%[0.5;1]
%[-0.5;2]
%[-1.5;-1.5]
%[-1;-1]
%[-1.5;0.5]
nav_p.xd = [1;0];%[1;0]; % Goal position / Desired position
nav_p.p = 2; % Generate obstacle set off p-norm
nav_p.rad_from_goal = 0.01;
nav_p.r1 = 0.5;
nav_p.r2 = 0.6;
nav_p.alpha = 0.2;
nav_p.beta = 1e-3;
nav_p.ctrl_num = 1;
nav_p.dim = 2;
nav_p.A = [0,1;1-2*nav_p.xd(1),-0.1];
nav_p.B = [0;1];
nav_p.Q = eye(nav_p.dim);
nav_p.R = 1;
N = 5e3; %timesteps
deltaT = 1e-2;
ctrl_bound = 4;
ctrl_multiplier = 1;
theta = 0*pi/180; % [Radian] CounterClockwise | Rotate after
stretch = [1;1];
epsilon = 0.001;


% Obstacle centers
num_obs = 1;
nav_p.c1 = [0;0];%[0; 0]; 
nav_p.c2 = [0; -2];
nav_p.c3 = [4; 0];


%% Density Function Formulation
syms x [nav_p.dim,1] real

[A1, A_inv] = transformationMatrix(theta, stretch, 2);

[K,P1,e] = lqr(nav_p.A,nav_p.B,nav_p.Q,nav_p.R);


% bump = formPNormBump(nav_p.r1,nav_p.r2, nav_p.c1, x, nav_p.p, true, A_inv)*...
%     formPNormBump(nav_p.r1,nav_p.r2, nav_p.c2, x, nav_p.p, true, A_inv)*...
%     formPNormBump(nav_p.r1,nav_p.r2, nav_p.c3, x, nav_p.p, true, A_inv);

bump = formFastInvBump(x, nav_p, num_obs);

g = 1/((x-nav_p.xd)'*P1*(x-nav_p.xd))^(nav_p.alpha);
density = g*bump;

optimize = false;

drift = [x2;...
         x1 - x1^3 - 0.1*x2];
matlabFunction(drift, 'File', 'functions/drift_f', 'Vars', {x}, 'Optimize', optimize);

matlabFunction(density, 'File', 'functions/density_f', 'Vars', {x}, 'Optimize', optimize);

div0 = divergence(drift*density,x);
matlabFunction(div0, 'File', 'functions/div0_f', 'Vars', {x}, 'Optimize', optimize);
div1 = divergence(nav_p.B(:,1)*density,x);
matlabFunction(div1, 'File', 'functions/div1_f', 'Vars', {x}, 'Optimize', optimize);




%% Linear programming


H_cdf = 2*eye(3*nav_p.ctrl_num+1);
f_cdf = zeros(3*nav_p.ctrl_num+1,1);
u_cdf = zeros(3*nav_p.ctrl_num+1,N);
x_cdf = zeros(nav_p.dim,N+1);
x_cdf(:,1) = nav_p.x0;
x_unctrl_cdf = zeros(nav_p.dim,N+1);
x_unctrl_cdf(:,1) = nav_p.x0;
x1_value = zeros(nav_p.dim,1);
x2_value = zeros(nav_p.dim,1);
u_cost = zeros(1,N);


% options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
% options = optimoptions('linprog','Algorithm','interior-point','Display','off');
options = optimoptions('quadprog','Display','off');

for i = 1:N
    if mod(i,100000) == 0
        i
    end
    
   
    if(norm(x_cdf(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
%         u_cdf(:,i) =  0;
        u_cdf(:,i)=-K*(x_cdf(:,i)-nav_p.xd);
%         u_cdf(:,i)=0;
%         break;
    else
        x1_value = x_cdf(:,i) + [epsilon;0];
        x2_value = x_cdf(:,i) + [0;epsilon];    

        A_cdf = [-div1_f(x_cdf(:,i)),0,0,0;...
            0,-div1_f(x1_value),0,0;...
            0,0,-div1_f(x2_value),0;...
            (-sum(nav_p.B(:,1)))/epsilon,nav_p.B(1,1)/epsilon,nav_p.B(2,1)/epsilon,1;...
            (sum(nav_p.B(:,1)))/epsilon,-nav_p.B(1,1)/epsilon,-nav_p.B(2,1)/epsilon,1];

        b_cdf = [div0_f(x_cdf(:,i))- nav_p.beta*density_f(x_cdf(:,i));...
                 div0_f(x1_value)- nav_p.beta*density_f(x1_value);...
                 div0_f(x2_value)- nav_p.beta*density_f(x2_value);...
                 nav_p.beta;...
                 nav_p.beta];
        lb_cdf = [-ctrl_bound*ones(3*nav_p.ctrl_num,1);1e-10];
        ub_cdf = [ctrl_bound*ones(3*nav_p.ctrl_num,1);Inf];
%         f_cdf = 2*K*(x_cdf(:,i)-nav_p.xd);


        u_cdf(:,i) = ctrl_multiplier*quadprog(H_cdf,f_cdf,A_cdf,b_cdf,[],[],lb_cdf,ub_cdf,[],options);
%         u_cdf(:,i) = ctrl_multiplier*linprog(f_cdf,A_cdf,b_cdf,[],[],lb_cdf,ub_cdf,options);

        
    end
    
    u_cost(:,i) = (u_cdf(1,i))^2;
    x_cdf(:,i+1) = x_cdf(:,i) + (drift_f(x_cdf(:,i)) + nav_p.B*u_cdf(1,i))*deltaT;
    x_unctrl_cdf(:,i+1) = x_unctrl_cdf(:,i) + drift_f(x_unctrl_cdf(:,i))*deltaT;
end

%% Plots (sriram)
colors = colororder;
blue = colors(1,:);
red = colors(2,:);
yellow = colors(3,:);
green = colors(5,:);
obsColor = [.7 .7 .7]; % Obstacle color -> Grey

figure(1)
hold on
%------------------------------ plot drift trajectory ------------------------------------
subplot(2,1,1)
Dom = [-2,2];
[X,Y] = meshgrid(Dom(1):0.25:Dom(2),Dom(1):0.25:Dom(2));
u = Y;
v = X - X.^3 - 0.1.*Y;
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color',obsColor);

%------------------------------ plot state trajectory ------------------------------------
state_traj_cdf = plot(x_cdf(1,:),x_cdf(2,:),'-', 'LineWidth', 2); hold on;
plot(nav_p.x0(1), nav_p.x0(2), 'o', 'MarkerSize',10, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
plot(nav_p.xd(1), nav_p.xd(2), 'o', 'MarkerSize',10, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;

% plot obstacles
gamma = (0:100-1)*(2*pi/100);
points = nav_p.c1 + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;

%plot settings
axes1 = gca;
box(axes1,'on');
axis(axes1,'square');
set(axes1,'FontSize',15,'LineWidth',2);
xlim([-2,2]);
ylim([-2,2]);
xlabel('$x_1\;(m)$','interpreter','latex', 'FontSize', 20);
ylabel('$x_2\;(m)$','interpreter','latex', 'FontSize', 20);

% lgd = legend('CDF','Uncontrolled', 'Start','Goal','Obstacles');
% lgd.FontSize = 14;
% set(lgd,'interpreter','latex')
% %title('Euclidean Space Trajectory');
% xlabel('$x_1$','interpreter','latex', 'FontSize', 20);
% ylabel('$x_2$','interpreter','latex', 'FontSize', 20);

%------------------------------ plot control values ------------------------------------
hold on
subplot(2,1,2)
t = 1:size(u_cdf,2);
t1 = 1:size(x_cdf,2);

plot(t*deltaT,u_cdf(1,:), 'LineWidth', 2); hold on;

%plot settings
axes2 = gca;
box(axes2,'on');
set(axes2,'FontSize',15,'LineWidth',2);
xlim([0,N*deltaT]);ylim([-4,4]);
xlabel("$\textrm{time},\;t\;(s)$",'interpreter','latex', 'FontSize', 20)
ylabel("$\textrm{control},u$",'interpreter','latex', 'FontSize', 20)
ctrl_lgd = legend('$u$','interpreter','latex');
ctrl_lgd.FontSize = 15;
xlim([0,10])

%------------------------------ other plots ------------------------------------
% figure()
% plot(t*deltaT,u_cdf(2,:), 'LineWidth', 2); hold on;
% hold off;
% title("u2 Control trajectory")
% xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
% ylabel("$u_2$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend('$u_2\; cdf$','interpreter','latex');
% ctrl_lgd.FontSize = 14;
% xlim([0,N*deltaT]);
% grid on;

% figure()
% plot(t,u_cdf(3,:), 'LineWidth', 2); hold on;
% hold off;
% title("u3 Control trajectory")
% xlabel("$Iteration$",'interpreter','latex', 'FontSize', 20)
% ylabel("$u_3$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend('$u_3\; cdf$','interpreter','latex');
% ctrl_lgd.FontSize = 14;
% xlim([0,N]);
% grid on;
% 
% figure()
% plot(t,u_cdf(4,:), 'LineWidth', 2); hold on;
% hold off;
% title("u4 Control trajectory")
% xlabel("$Iteration$",'interpreter','latex', 'FontSize', 20)
% ylabel("$u_4$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend('$u_4\; cdf$','interpreter','latex');
% ctrl_lgd.FontSize = 14;
% xlim([0,N]);
% grid on;
% 
% figure()
% plot(t,u_cdf(5,:), 'LineWidth', 2); hold on;
% hold off;
% title("u5 Control trajectory")
% xlabel("$Iteration$",'interpreter','latex', 'FontSize', 20)
% ylabel("$u_5$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend('$u_5\; cdf$','interpreter','latex');
% ctrl_lgd.FontSize = 14;
% xlim([0,N]);
% grid on;
% 
% figure()
% plot(t,u_cdf(6,:), 'LineWidth', 2); hold on;
% hold off;
% title("u6 Control trajectory")
% xlabel("$Iteration$",'interpreter','latex', 'FontSize', 20)
% ylabel("$u_6$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend('$u_6\; cdf$','interpreter','latex');
% ctrl_lgd.FontSize = 14;
% xlim([0,N]);
% grid on;

% figure()
% plot(t,u_cost(1,:), 'LineWidth', 2); hold on;
% hold off;
% title("Control Cost")
% xlabel("$Iteration$",'interpreter','latex', 'FontSize', 20)
% ylabel("$u\;cost$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd.FontSize = 14;
% xlim([0,N]);
% grid on;

