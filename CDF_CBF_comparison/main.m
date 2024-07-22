%% Obstacle Avoidance control using quadratic program
clc; clear; close all; 
mkdir('functions');
addpath('./functions');

%% Initialization
nav_p.x0 = [-4.9;0.1]; % Start position
nav_p.xd = [4.9;0]; % Goal position / Desired position
nav_p.p = 2; % Generate obstacle set off p-norm
nav_p.rad_from_goal = 0.1;
nav_p.r1 = 1;
nav_p.r2_1 = 2;
nav_p.r2_2 = 3;
nav_p.r2_3 = 4;
nav_p.alpha = 1;
nav_p.beta = 1e-3;

nav_p.gamma1 = 0.7;  %CBF exponential coefficient 1
nav_p.gamma2 = 0.5;  %CBF exponential coefficient 2
nav_p.gamma3 = 0.3;  %CBF exponential coefficient 3
nav_p.lambda = 0.1; % CLF exponential coefficient

nav_p.ctrl_num = 2;
nav_p.dim = 2;
nav_p.A = zeros(2);%[0 0;-1 0];
nav_p.B = eye(nav_p.dim);
nav_p.Q = eye(nav_p.dim);
nav_p.R = 1;
N = 1.1e3; %timesteps
deltaT = 1e-2;
ctrl_bound = 2;
ctrl_multiplier = 1;
theta = 0*pi/180; % [Radian] CounterClockwise | Rotate after
stretch = [1;1];
epsilon = 0.001;


% Obstacle centers
num_obs = 3;
nav_p.c1 = [0;0]; 



%% Density Function Formulation
syms x [nav_p.dim,1] real

[A1, A_inv] = transformationMatrix(theta, stretch, 2);

[K,P1,e] = lqr(nav_p.A,nav_p.B,nav_p.Q,nav_p.R);

h1 = norm(x-nav_p.c1) - nav_p.r1;
% h2 = norm(x-nav_p.c2) - nav_p.r1;
% h3 = norm(x-nav_p.c3) - nav_p.r1;

grad_h1 = gradient(h1, x);
% grad_h2 = gradient(h2, x);
% grad_h3 = gradient(h3, x);

V = (x-nav_p.xd)'*P1*(x-nav_p.xd);

grad_V = gradient(V, x);

optimize = false;

matlabFunction(h1, 'File', 'functions/h1_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(grad_h1, 'File', 'functions/grad_h1_f', 'Vars', {x}, 'Optimize', optimize);
% matlabFunction(h2, 'File', 'functions/h2_f', 'Vars', {x}, 'Optimize', optimize);
% matlabFunction(grad_h2, 'File', 'functions/grad_h2_f', 'Vars', {x}, 'Optimize', optimize);
% matlabFunction(h3, 'File', 'functions/h3_f', 'Vars', {x}, 'Optimize', optimize);
% matlabFunction(grad_h3, 'File', 'functions/grad_h3_f', 'Vars', {x}, 'Optimize', optimize);

matlabFunction(V, 'File', 'functions/V_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(grad_V, 'File', 'functions/grad_V_f', 'Vars', {x}, 'Optimize', optimize);

bump1 = formPNormBump(nav_p.r1,nav_p.r2_1, nav_p.c1, x, nav_p.p, true, A_inv);
bump2 = formPNormBump(nav_p.r1,nav_p.r2_2, nav_p.c1, x, nav_p.p, true, A_inv);
bump3 = formPNormBump(nav_p.r1,nav_p.r2_3, nav_p.c1, x, nav_p.p, true, A_inv);


% bump = formFastInvBump(x, nav_p, num_obs);

g = 1/((x-nav_p.xd)'*P1*(x-nav_p.xd))^(nav_p.alpha);
density1 = g*bump1;
density2 = g*bump2;
density3 = g*bump3;


matlabFunction(density1, 'File', 'functions/density1_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(density2, 'File', 'functions/density2_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(density3, 'File', 'functions/density3_f', 'Vars', {x}, 'Optimize', optimize);


div0_1 = divergence(nav_p.A*x*density1,x);
matlabFunction(div0_1, 'File', 'functions/div0_1_f', 'Vars', {x}, 'Optimize', optimize);
div1_1 = divergence(nav_p.B(:,1)*density1,x);
matlabFunction(div1_1, 'File', 'functions/div1_1_f', 'Vars', {x}, 'Optimize', optimize);
div2_1 = divergence(nav_p.B(:,2)*density1,x);
matlabFunction(div2_1, 'File', 'functions/div2_1_f', 'Vars', {x}, 'Optimize', optimize);

div0_2 = divergence(nav_p.A*x*density2,x);
matlabFunction(div0_2, 'File', 'functions/div0_2_f', 'Vars', {x}, 'Optimize', optimize);
div1_2 = divergence(nav_p.B(:,1)*density2,x);
matlabFunction(div1_2, 'File', 'functions/div1_2_f', 'Vars', {x}, 'Optimize', optimize);
div2_2 = divergence(nav_p.B(:,2)*density2,x);
matlabFunction(div2_2, 'File', 'functions/div2_2_f', 'Vars', {x}, 'Optimize', optimize);

div0_3 = divergence(nav_p.A*x*density3,x);
matlabFunction(div0_3, 'File', 'functions/div0_3_f', 'Vars', {x}, 'Optimize', optimize);
div1_3 = divergence(nav_p.B(:,1)*density3,x);
matlabFunction(div1_3, 'File', 'functions/div1_3_f', 'Vars', {x}, 'Optimize', optimize);
div2_3 = divergence(nav_p.B(:,2)*density3,x);
matlabFunction(div2_3, 'File', 'functions/div2_3_f', 'Vars', {x}, 'Optimize', optimize);



%% Linear programming

no_decision_variable = 2*nav_p.ctrl_num+1;

H_cdf = 2*blkdiag(1,1,0,0,1);
f_cdf = zeros(no_decision_variable,1);

u_cdf_1 = zeros(no_decision_variable,N);
u_cdf_2 = zeros(no_decision_variable,N);
u_cdf_3 = zeros(no_decision_variable,N);

x_cdf_1 = zeros(nav_p.dim,N+1);
x_cdf_1(:,1) = nav_p.x0;
x_cdf_2 = zeros(nav_p.dim,N+1);
x_cdf_2(:,1) = nav_p.x0;
x_cdf_3 = zeros(nav_p.dim,N+1);
x_cdf_3(:,1) = nav_p.x0;

x1_value = zeros(nav_p.dim,1);
x2_value = zeros(nav_p.dim,1);
% u_cost = zeros(1,N);

H_cbf = 2*eye(nav_p.ctrl_num);
f_cbf = zeros(nav_p.ctrl_num,1);

u_cbf_1 = zeros(nav_p.ctrl_num,N);
u_cbf_2 = zeros(nav_p.ctrl_num,N);
u_cbf_3 = zeros(nav_p.ctrl_num,N);

x_cbf_1 = zeros(nav_p.dim,N+1);
x_cbf_1(:,1) = nav_p.x0;
x_cbf_2 = zeros(nav_p.dim,N+1);
x_cbf_2(:,1) = nav_p.x0;
x_cbf_3 = zeros(nav_p.dim,N+1);
x_cbf_3(:,1) = nav_p.x0;
% u_cost = zeros(1,N);


% options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
% options = optimoptions('linprog','Algorithm','interior-point','Display','off');
options = optimoptions('quadprog','Display','off');

for i = 1:N
    
    
    if mod(i,1e4) == 0
        i
    end
    
    if(norm(x_cbf_1(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
        u_cbf_1(:,i) =  0; 
    else
    
        %     Three Obstacle with Lyapunov function

        A_cbf_1 = [grad_V_f(x_cbf_1(:,i))'*nav_p.B(:,1),grad_V_f(x_cbf_1(:,i))'*nav_p.B(:,2);...
            -grad_h1_f(x_cbf_1(:,i))'*nav_p.B(:,1),-grad_h1_f(x_cbf_1(:,i))'*nav_p.B(:,2)];

        b_cbf_1 = [-grad_V_f(x_cbf_1(:,i))'*nav_p.A*x_cbf_1(:,i) - nav_p.lambda*V_f(x_cbf_1(:,i));...
                grad_h1_f(x_cbf_1(:,i))'*nav_p.A*x_cbf_1(:,i) + nav_p.gamma1*h1_f(x_cbf_1(:,i))];

        lb_cbf = -ctrl_bound*ones(nav_p.ctrl_num,1);
        ub_cbf = ctrl_bound*ones(nav_p.ctrl_num,1);
        f_cbf = 2*K*(x_cbf_1(:,i)-nav_p.xd);

        u_cbf_1(:,i) = ctrl_multiplier*quadprog(H_cbf,f_cbf,A_cbf_1,b_cbf_1,[],[],lb_cbf,ub_cbf,[],options);
    end
    if(norm(x_cbf_2(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
        u_cbf_2(:,i) =  0; 
    else
    
        %     Three Obstacle with Lyapunov function

        A_cbf_2 = [grad_V_f(x_cbf_2(:,i))'*nav_p.B(:,1),grad_V_f(x_cbf_2(:,i))'*nav_p.B(:,2);...
            -grad_h1_f(x_cbf_2(:,i))'*nav_p.B(:,1),-grad_h1_f(x_cbf_2(:,i))'*nav_p.B(:,2)];

        b_cbf_2 = [-grad_V_f(x_cbf_2(:,i))'*nav_p.A*x_cbf_2(:,i) - nav_p.lambda*V_f(x_cbf_2(:,i));...
                grad_h1_f(x_cbf_2(:,i))'*nav_p.A*x_cbf_2(:,i) + nav_p.gamma2*h1_f(x_cbf_2(:,i))];

        lb_cbf = -ctrl_bound*ones(nav_p.ctrl_num,1);
        ub_cbf = ctrl_bound*ones(nav_p.ctrl_num,1);
        f_cbf = 2*K*(x_cbf_2(:,i)-nav_p.xd);

        u_cbf_2(:,i) = ctrl_multiplier*quadprog(H_cbf,f_cbf,A_cbf_2,b_cbf_2,[],[],lb_cbf,ub_cbf,[],options);
    end
    if(norm(x_cbf_3(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
        u_cbf_3(:,i) =  0; 
    else
    
        %     Three Obstacle with Lyapunov function

        A_cbf_3 = [grad_V_f(x_cbf_3(:,i))'*nav_p.B(:,1),grad_V_f(x_cbf_3(:,i))'*nav_p.B(:,2);...
            -grad_h1_f(x_cbf_3(:,i))'*nav_p.B(:,1),-grad_h1_f(x_cbf_3(:,i))'*nav_p.B(:,2)];

        b_cbf_3 = [-grad_V_f(x_cbf_3(:,i))'*nav_p.A*x_cbf_3(:,i) - nav_p.lambda*V_f(x_cbf_3(:,i));...
                grad_h1_f(x_cbf_3(:,i))'*nav_p.A*x_cbf_3(:,i) + nav_p.gamma3*h1_f(x_cbf_3(:,i))];

        lb_cbf = -ctrl_bound*ones(nav_p.ctrl_num,1);
        ub_cbf = ctrl_bound*ones(nav_p.ctrl_num,1);
        f_cbf = 2*K*(x_cbf_3(:,i)-nav_p.xd);

        u_cbf_3(:,i) = ctrl_multiplier*quadprog(H_cbf,f_cbf,A_cbf_3,b_cbf_3,[],[],lb_cbf,ub_cbf,[],options);
    end
    
    
    
    
    
    if(norm(x_cdf_1(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
        u_cdf_1(1:2,i) =  -K*(x_cdf_1(:,i)-nav_p.xd); 
    else
        x1_value = x_cdf_1(:,i) + nav_p.B(:,1)*deltaT*density1_f(x_cdf_1(:,i));
        x2_value = x_cdf_1(:,i) + nav_p.B(:,1)*deltaT*density1_f(x_cdf_1(:,i));    

        A_cdf_1 = [-div1_1_f(x_cdf_1(:,i)),-div2_1_f(x_cdf_1(:,i)),0,0,1;...
                 0,0,-div1_1_f(x1_value),0,1;...
                 0,0,0,-div2_1_f(x2_value),1;...
                 1/deltaT,1/deltaT,-1/deltaT,-1/deltaT,0];

        b_cdf_1 = [div0_1_f(x_cdf_1(:,i));div0_1_f(x1_value);div0_1_f(x2_value);0];
        lb_cdf = [-ctrl_bound*ones(1*nav_p.ctrl_num,1);-1e5;-1e5;1e-5];
        ub_cdf = [ctrl_bound*ones(1*nav_p.ctrl_num,1);1e5;1e5;1e5];
        
        f_cdf = [2*K*(x_cdf_1(:,i)-nav_p.xd);0;0;0];


        u_cdf_1(:,i) = quadprog(H_cdf,f_cdf,A_cdf_1,b_cdf_1,[],[],lb_cdf,ub_cdf,[],options);
        
    end
    if(norm(x_cdf_2(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
        u_cdf_2(1:2,i) =  -K*(x_cdf_2(:,i)-nav_p.xd); 
    else
        x1_value = x_cdf_2(:,i) + nav_p.B(:,1)*deltaT*density2_f(x_cdf_2(:,i));
        x2_value = x_cdf_2(:,i) + nav_p.B(:,1)*deltaT*density2_f(x_cdf_2(:,i));    

        A_cdf_2 = [-div1_2_f(x_cdf_2(:,i)),-div2_2_f(x_cdf_2(:,i)),0,0,1;...
                 0,0,-div1_2_f(x1_value),0,1;...
                 0,0,0,-div2_2_f(x2_value),1;...
                 1/deltaT,1/deltaT,-1/deltaT,-1/deltaT,0];

        b_cdf_2 = [div0_2_f(x_cdf_2(:,i));div0_2_f(x1_value);div0_2_f(x2_value);0];
        lb_cdf = [-ctrl_bound*ones(1*nav_p.ctrl_num,1);-1e5;-1e5;1e-5];
        ub_cdf = [ctrl_bound*ones(1*nav_p.ctrl_num,1);1e5;1e5;1e5];
        
        f_cdf = [2*K*(x_cdf_2(:,i)-nav_p.xd);0;0;0];


        u_cdf_2(:,i) = quadprog(H_cdf,f_cdf,A_cdf_2,b_cdf_2,[],[],lb_cdf,ub_cdf,[],options);
        
    end
    if(norm(x_cdf_3(:,i)-nav_p.xd)<nav_p.rad_from_goal)

        % LQR Feedback Gain
        u_cdf_3(1:2,i) =  -K*(x_cdf_3(:,i)-nav_p.xd); 
    else
        x1_value = x_cdf_3(:,i) + nav_p.B(:,1)*deltaT*density3_f(x_cdf_3(:,i));
        x2_value = x_cdf_3(:,i) + nav_p.B(:,1)*deltaT*density3_f(x_cdf_3(:,i));    

        A_cdf_3 = [-div1_3_f(x_cdf_3(:,i)),-div2_3_f(x_cdf_3(:,i)),0,0,1;...
                 0,0,-div1_3_f(x1_value),0,1;...
                 0,0,0,-div2_3_f(x2_value),1;...
                 1/deltaT,1/deltaT,-1/deltaT,-1/deltaT,0];

        b_cdf_3 = [div0_3_f(x_cdf_3(:,i));div0_3_f(x1_value);div0_3_f(x2_value);0];
        lb_cdf = [-ctrl_bound*ones(1*nav_p.ctrl_num,1);-1e5;-1e5;1e-5];
        ub_cdf = [ctrl_bound*ones(1*nav_p.ctrl_num,1);1e5;1e5;1e5];
        
        f_cdf = [2*K*(x_cdf_3(:,i)-nav_p.xd);0;0;0];


        u_cdf_3(:,i) = quadprog(H_cdf,f_cdf,A_cdf_3,b_cdf_3,[],[],lb_cdf,ub_cdf,[],options);
        
    end
    
    
    
    x_cdf_1(:,i+1) = x_cdf_1(:,i) + (nav_p.A*x_cdf_1(:,i) + nav_p.B*u_cdf_1(1:2,i))*deltaT;
    x_cbf_1(:,i+1) = x_cbf_1(:,i) + (nav_p.A*x_cbf_1(:,i) + nav_p.B*u_cbf_1(1:2,i))*deltaT;
    x_cdf_2(:,i+1) = x_cdf_2(:,i) + (nav_p.A*x_cdf_2(:,i) + nav_p.B*u_cdf_2(1:2,i))*deltaT;
    x_cbf_2(:,i+1) = x_cbf_2(:,i) + (nav_p.A*x_cbf_2(:,i) + nav_p.B*u_cbf_2(1:2,i))*deltaT;
    x_cdf_3(:,i+1) = x_cdf_3(:,i) + (nav_p.A*x_cdf_3(:,i) + nav_p.B*u_cdf_3(1:2,i))*deltaT;
    x_cbf_3(:,i+1) = x_cbf_3(:,i) + (nav_p.A*x_cbf_3(:,i) + nav_p.B*u_cbf_3(1:2,i))*deltaT;
end

%% Plot
close all;
obsColor = [.7 .7 .7]; % Obstacle color -> Grey

figure()
plot(x_cbf_1(1,:),x_cbf_1(2,:),'r--', 'LineWidth', 1.5); hold on;
plot(x_cbf_2(1,:),x_cbf_2(2,:),'Color','#A0F','LineStyle','--', 'LineWidth', 1.5); hold on;
plot(x_cbf_3(1,:),x_cbf_3(2,:),'b--', 'LineWidth', 1.5); hold on;
plot(x_cdf_1(1,:),x_cdf_1(2,:),'r', 'LineWidth', 1.5); hold on;
plot(x_cdf_2(1,:),x_cdf_2(2,:),'Color','#A0F','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(x_cdf_3(1,:),x_cdf_3(2,:),'b', 'LineWidth', 1.5); hold on;
plot(nav_p.x0(1), nav_p.x0(2), 'ok', 'MarkerSize',7, 'MarkerFaceColor','black'); hold on;
plot(nav_p.xd(1), nav_p.xd(2), 'og', 'MarkerSize',7, 'MarkerFaceColor','green'); hold on;
dummy_marker = plot(NaN,NaN, 'o','MarkerSize', 10, 'MarkerEdgeColor',...
            'black', 'MarkerFaceColor',obsColor, 'LineWidth', 1.5); % For legend as rectangular object can't be defined as a legend

gamma = (0:100-1)*(2*pi/100);
points = nav_p.c1 + A1*[nav_p.r2_3*cos(gamma);nav_p.r2_3*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', [0.3010 0.7450 0.9330], 'LineWidth', 1, 'EdgeColor', 'none','FaceAlpha', 0.2); hold on;

points = nav_p.c1 + A1*[nav_p.r2_2*cos(gamma);nav_p.r2_2*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', [ 0.4940 0.1840 0.5560], 'LineWidth', 1, 'EdgeColor', 'none','FaceAlpha', 0.2); hold on;

points = nav_p.c1 + A1*[nav_p.r2_1*cos(gamma);nav_p.r2_1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1, 'EdgeColor', 'none','FaceAlpha', 0.2); hold on;


points = nav_p.c1 + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 1.5, 'FaceAlpha', 1.0); hold on;
grid off;
hold off;
xlim([-5,5]);
ylim([-5,5]);
axis square

lgd = legend('CBF, $\gamma = 0.7$','CBF, $\gamma = 0.5$', 'CBF, $\gamma = 0.3$','CDF, $s=2$','CDF, $s=3$','CDF, $s=4$');
lgd.FontSize = 20;
lgd.NumColumns = 2;
set(lgd,'interpreter','latex')
%title('Euclidean Space Trajectory');
xlabel('position $x_1$ [m]','interpreter','latex', 'FontSize', 30);
ylabel('position $x_2$ [m]','interpreter','latex', 'FontSize', 30);


t = 1:size(u_cdf_1,2);
t1 = 1:size(x_cdf_1,2);

figure()
subplot(3,1,1)
plot(t*deltaT,u_cdf_1(1,:),'Color','#0072BD','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cdf_1(2,:),'Color','#D95319','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cbf_1(1,:),'Color','#D95319','LineStyle','--', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cbf_1(2,:),'Color','#0072BD','LineStyle','--', 'LineWidth', 1.5); hold on;

hold off;
title("Control plot for CBF with $\gamma=2$ and CDF with $s=2$",'interpreter','latex', 'FontSize', 45)
% xlabel("$time\;[s]$",'interpreter','latex', 'FontSize', 20)
ylabel("$\textbf{u}$",'interpreter','latex', 'FontSize', 45)
ctrl_lgd = legend('CDF $u_1$','CDF $u_2$','CBF $u_1$','CBF $u_2$','interpreter','latex');
ctrl_lgd.FontSize = 20;
ctrl_lgd.NumColumns = 2;
xlim([0,N*deltaT]);
ylim([-2.2,2.2]);
grid off;

subplot(3,1,2)
plot(t*deltaT,u_cdf_2(1,:),'Color','#0072BD','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cdf_2(2,:),'Color','#D95319','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cbf_2(1,:),'Color','#D95319','LineStyle','--', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cbf_2(2,:),'Color','#0072BD','LineStyle','--', 'LineWidth', 1.5); hold on;

hold off;
title("Control plot for CBF with $\gamma=0.5$ and CDF with $s=3$",'interpreter','latex', 'FontSize', 45)
% xlabel("$time\;[s]$",'interpreter','latex', 'FontSize', 20)
ylabel("$\textbf{u}$",'interpreter','latex', 'FontSize', 45)
ctrl_lgd = legend('CDF $u_1$','CDF $u_2$','CBF $u_1$','CBF $u_2$','interpreter','latex');
ctrl_lgd.FontSize = 20;
ctrl_lgd.NumColumns = 2;
xlim([0,N*deltaT]);
ylim([-2.2,2.2]);
grid off;

subplot(3,1,3)
plot(t*deltaT,u_cdf_3(1,:),'Color','#0072BD','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cdf_3(2,:),'Color','#D95319','LineStyle','-', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cbf_3(1,:),'Color','#D95319','LineStyle','--', 'LineWidth', 1.5); hold on;
plot(t*deltaT,u_cbf_3(2,:),'Color','#0072BD','LineStyle','--', 'LineWidth', 1.5); hold on;

hold off;
title("Control plot for CBF with $\gamma=0.2$ and CDF with $s=4$",'interpreter','latex', 'FontSize', 45)
xlabel("$time\;[s]$",'interpreter','latex', 'FontSize', 45)
ylabel("$\textbf{u}$",'interpreter','latex', 'FontSize', 45)
ctrl_lgd = legend('CDF $u_1$','CDF $u_2$','CBF $u_1$','CBF $u_2$','interpreter','latex');
ctrl_lgd.FontSize = 20;
ctrl_lgd.NumColumns = 2;
xlim([0,N*deltaT]);
ylim([-2.2,2.2]);
grid off;




