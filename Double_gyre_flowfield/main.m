%% Obstacle Avoidance control using quadratic program
clc; clear; close all; 
mkdir('functions');
addpath('./functions');

%% Initialization

% start at one eqb and end at the other eqb point
% note that [0.5,0.5] does not work well as it needs very fine discretization
nav_p.x0 = [1.5;0.5]; % Start position
nav_p.xd = [0.51;0.51]; % Goal position / Desired position
deltaT = 1e-4; % needs very small time steps to work well
deltaT1 = 1e-4;
ctrl_bound = 2*pi; % play with control bounds to get differnt solutions


% other params setup
nav_p.p = 2; 
nav_p.rad_from_goal = 0.01;
nav_p.r1 = 0.25;
nav_p.r2 = 0.3;
nav_p.alpha = 0.2;
nav_p.beta = 1e-3;
nav_p.ctrl_num = 2;
nav_p.dim = 2;
nav_p.A = [-pi 0; 0 pi];
nav_p.B = 2*eye(2);
nav_p.Q = eye(nav_p.dim);
nav_p.R = 1;
N = 4.2e4; 

theta = 0*pi/180; 
stretch = [1;1];
epsilon = 0.001;

% Obstacle centers
num_obs = 1;
nav_p.c1 = [1; 0]; 

%% Density Function Formulation
syms x [nav_p.dim,1] real

[A1, A_inv] = transformationMatrix(theta, stretch, 2);
[K,P1,e] = lqr(nav_p.A,nav_p.B,nav_p.Q,nav_p.R);
optimize = false;
hamiltonian = -sin(pi*x1)*sin(pi*x2);
matlabFunction(hamiltonian, 'File', 'functions/hamiltonian_f', 'Vars', {x}, 'Optimize', optimize);
bump = formFastInvBump(x, nav_p, num_obs);
g = 1/((x-nav_p.xd)'*eye(2)*(x-nav_p.xd))^(nav_p.alpha);
% g = 1/(hamiltonian)^(2);
density = g*bump;

drift = [-pi*sin(pi*x1)*cos(pi*x2);...
         pi*sin(pi*x2)*cos(pi*x1)];

matlabFunction(drift, 'File', 'functions/drift_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(density, 'File', 'functions/density_f', 'Vars', {x}, 'Optimize', optimize);

div0 = divergence(drift*density,x);  %divergence(f(x)*\rho)
matlabFunction(div0, 'File', 'functions/div0_f', 'Vars', {x}, 'Optimize', optimize);
div1 = divergence(nav_p.B(:,1)*density,x); %divergence(g1(x)*\rho)
matlabFunction(div1, 'File', 'functions/div1_f', 'Vars', {x}, 'Optimize', optimize);
div2 = divergence(nav_p.B(:,2)*density,x); %divergence(g2(x)*\rho)
matlabFunction(div2, 'File', 'functions/div2_f', 'Vars', {x}, 'Optimize', optimize);

%% Linear programming

% add is contrl multiplier gain, fixes most of the issues
gain = 2; 

no_decision_variable = 2*nav_p.ctrl_num+1;
H_cdf = 2*blkdiag(1,1,0,0,1);
f_cdf = zeros(no_decision_variable,1);
u_cdf = NaN*ones(no_decision_variable,N);
x_cdf = NaN*ones(nav_p.dim,N+1);
x_cdf(:,1) = nav_p.x0;
x1_value = zeros(nav_p.dim,1);
x2_value = zeros(nav_p.dim,1);
options = optimoptions('quadprog','Display','off','MaxIterations',1e6,'Algorithm','interior-point-convex');

for i = 1:N
    if mod(i,1e4) == 0
        i
    end
  
    if(norm(x_cdf(:,i)-nav_p.xd)<nav_p.rad_from_goal)
        % exit out of the simulation when you are close to LQR doamin
        disp('reached near target')
        u_cdf(1:2,i) =  [0;0]; 
        break 
    else
        % get control using QP formulation
        x1_value = x_cdf(:,i) + nav_p.B(:,1)*deltaT*density_f(x_cdf(:,i)); %x1 which is 'x' propogated along the flow of g1(x)*\rho 
        x2_value = x_cdf(:,i) + nav_p.B(:,2)*deltaT*density_f(x_cdf(:,i)); %x2 which is 'x' propogated along the flow of g2(x)*\rho 
        A_cdf = [-div1_f(x_cdf(:,i)),-div2_f(x_cdf(:,i)),0,0,1;...
                 0,0,-div1_f(x1_value),0,1;...
                 0,0,0,-div2_f(x2_value),1;...
                 1/deltaT1,1/deltaT1,-1/deltaT1,-1/deltaT1,0];           
        b_cdf = [div0_f(x_cdf(:,i));div0_f(x1_value);div0_f(x2_value);0];
        lb_cdf = [-ctrl_bound*ones(1*nav_p.ctrl_num,1);-1e5;-1e5;1e-5];
        ub_cdf = [ctrl_bound*ones(1*nav_p.ctrl_num,1);1e5;1e5;1e5];
        
        % add a control multiplier gain to get better control
        u_cdf(:,i) = gain*quadprog(H_cdf,f_cdf,A_cdf,b_cdf,[],[],lb_cdf,ub_cdf,[],options);  
    end
    
    % euler update
    x_cdf(:,i+1) = x_cdf(:,i) + (drift_f(x_cdf(:,i)) + nav_p.B*u_cdf(1:2,i))*deltaT;
end
%%
% for i =1:N
%     x_cdf(:,i+1) = x_cdf(:,i) + (drift_f(x_cdf(:,i)) + nav_p.B*gain*u_cdf(1:2,i))*deltaT;
% end

%% Plot
close all;
colors = colororder;
blue = colors(1,:);
red = colors(2,:);
yellow = colors(3,:);
green = colors(5,:);
obsColor = [.7 .7 .7]; % Obstacle color -> Grey

figure()
% plot the drift
Dom_x = [0 2];
Dom_y = [0 1];
[X,Y] = meshgrid(Dom_x(1):0.25:Dom_x(2),Dom_y(1):0.25:Dom_y(2));
u = -pi.*sin(pi.*X).*cos(pi.*Y);
v =  pi.*sin(pi.*Y).*cos(pi.*X);
l = streamslice(X,Y,u,v); hold on;
set(l,'LineWidth',1)
set(l,'Color',obsColor);
hold on

% plot start goal and traj
state_traj_cdf = plot(x_cdf(1,:),x_cdf(2,:),'Color',blue, 'LineWidth', 2); hold on;
plot(nav_p.x0(1), nav_p.x0(2), 'o', 'MarkerSize',15, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
plot(nav_p.xd(1), nav_p.xd(2), 'o', 'MarkerSize',15, 'MarkerFaceColor',green,'MarkerEdgeColor',green); hold on;

% plot obs
gamma = (0:100-1)*(2*pi/100);
points = nav_p.c1 + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;

% plot settings
grid off;
axes1 = gca;
box(axes1,'on');
set(axes1,'FontSize',15,'LineWidth',2);
xlim([0,2]);
% ylim([-0.5,1.5]); axis square
ylim([0,1]);
xlabel('position $x_1$ [m]','interpreter','latex', 'FontSize', 20);
ylabel('position $x_2$ [m]','interpreter','latex', 'FontSize', 20);

%% control plots
t = 1:size(u_cdf,2);
t1 = 1:size(x_cdf,2);
figure()
plot(t*deltaT,u_cdf(1,:),'Color',blue, 'LineWidth', 2); hold on;
plot(t*deltaT,u_cdf(2,:),'Color',red,'LineWidth', 2, 'LineStyle','--'); hold on;
hold off;
% title("u1 Control trajectory")
xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$\textbf{u}$",'interpreter','latex', 'FontSize', 20)
ctrl_lgd = legend('$u_1$','$u_2$','interpreter','latex');
ctrl_lgd.FontSize = 14;
xlim([0,N*deltaT]);
grid on;