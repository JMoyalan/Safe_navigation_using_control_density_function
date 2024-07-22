%% Obstacle Avoidance control using quadratic program
clc; clear; close all; 
mkdir('functions');
addpath('./functions');

%% Initialization
nav_p.x0 = [-0.5;0;0;0]; % Start position
nav_p.p = 2; % Generate obstacle set off p-norm
nav_p.rad_from_goal = 0.1;
nav_p.alpha = 30;  %30(30) for x0 at +ve side and safe(sensing) region || 30(30) for x0 at -ve side and safe(sensing) region
nav_p.beta = 1e-3;
nav_p.ctrl_num = 1;
nav_p.dim = 4;
nav_p.disp_max_h = 0.9;%0.9;
nav_p.disp_max_s = 0.7;%0.55;
nav_p.acc_max = 0.5*9.81;   % 0.5(0.5) for x0 at +ve side and safe(sensing) region || 0.3(0.3) for x0 at -ve side and safe(sensing) region
nav_p.M = 1589;
nav_p.a = 1.57;
nav_p.b = 1.05;
nav_p.Cf = 90000;
nav_p.Cr = 60000;
nav_p.Iz = 1765;
nav_p.c = [0;0];
nav_p.L = 20;
nav_p.r = 300;%207;%107;%50;
nav_p.ctrl_bound = Inf;
N =1000; %timesteps
kappa = 11;%0.6
rng(0);
random_iter = 10;
alpha = 0.5;
Cr_max = 0000;
Cr_max1 = 5000;
rand_Cr = (2*Cr_max1*rand(random_iter,1) - Cr_max1);


delta_bound = 100;  
% nav_p.Q = 1*eye(nav_p.dim);
nav_p.Q = [1,0,0,0;...
           0,1,0,0;...
           0,0,1,0;...
           0,0,0,1];
nav_p.R =1500;     %usually we need less R value for +ive value of initial lateral displacement like the value 1
deltaT = 1e-2;
% ctrl_bound = 5;
ctrl_multiplier = 1;
epsilon = 1e-2;

%%
v02 = ((-2)/(nav_p.M))*(((((nav_p.a^2)*nav_p.Cf + (nav_p.b^2)*(nav_p.Cr+Cr_max))*((nav_p.Cr+Cr_max)+nav_p.Cf))/((nav_p.b*(nav_p.Cr+Cr_max)) - (nav_p.a*nav_p.Cf))) - (nav_p.b*(nav_p.Cr+Cr_max)) + (nav_p.a*nav_p.Cf));
v0 = sqrt(v02);
x_2 = 0;
x_3 = (((nav_p.a^2)*nav_p.Cf + (nav_p.b^2)*(nav_p.Cr+Cr_max))*(v0/nav_p.r))/(((nav_p.a*nav_p.Cf)-(nav_p.b*(nav_p.Cr+Cr_max)))*v0);
x_4 = v0/nav_p.r;

nav_p.v0 = v0;%20;%15;%15;
rd = nav_p.v0/nav_p.r;
%%



nav_p.A = [0,1,0,-nav_p.L;...
           0,-(2*((nav_p.Cr+Cr_max)+nav_p.Cf))/(nav_p.M*nav_p.v0),(2*((nav_p.Cr+Cr_max)+nav_p.Cf))/(nav_p.M),((2*((nav_p.b*(nav_p.Cr+Cr_max)) - (nav_p.a*nav_p.Cf)))/(nav_p.M*nav_p.v0))-(2*nav_p.v0);...
           0,0,0,-1;...
           0,(2*((nav_p.b*(nav_p.Cr+Cr_max)) - (nav_p.a*nav_p.Cf)))/(nav_p.Iz*nav_p.v0),(-2*((nav_p.b*(nav_p.Cr+Cr_max)) - (nav_p.a*nav_p.Cf)))/(nav_p.Iz),(-2*((nav_p.a^2)*nav_p.Cf + (nav_p.b^2)*(nav_p.Cr+Cr_max)))/(nav_p.Iz*nav_p.v0)];
               
        
nav_p.B = [0;2*nav_p.Cf/nav_p.M;0;(2*nav_p.a*nav_p.Cf)/(nav_p.Iz)];

nav_p.C = [nav_p.L;nav_p.v0;1;0]*rd;


r = rank(ctrb(nav_p.A,nav_p.B));


% Obstacle centers
% num_obs = 1;
% nav_p.c1 = [-5; 0.1]; 
% nav_p.c2 = [0; -2];
% nav_p.c3 = [4; 0];

% nav_p.c1 = [4; 4; 4; 4]; 

% r = rank(ctrb(nav_p.A,nav_p.B))


%% Density Function Formulation
syms x [nav_p.dim,1] real


[K1,P1,e1] = lqr(nav_p.A,nav_p.B,nav_p.Q,nav_p.R);
K1 = 3*K1;

% nav_p.A = nav_p.A-nav_p.B*K; % For stabilising nav_p.A
% r = rank(ctrb(nav_p.A,nav_p.B))

% bump = formPNormBump(nav_p.r1,nav_p.r2, nav_p.c1, x, nav_p.p, true, A_inv)*...
%     formPNormBump(nav_p.r1,nav_p.r2, nav_p.c2, x, nav_p.p, true, A_inv)*...
%     formPNormBump(nav_p.r1,nav_p.r2, nav_p.c3, x, nav_p.p, true, A_inv);

% bump = formFastInvBump(x, nav_p, num_obs);
% bump = nav_p.disp_max_h - sign(x2)*x1 - 0.5*(x2^2)/(nav_p.acc_max);

bump = LaneKeepingBump(nav_p.disp_max_h,nav_p.disp_max_s,nav_p.acc_max, [x1;x2]);

% for i = 1:1e20
%     P2 = rand(4);
%     P2 = P2+P2';
%     if not(any(eig(P2)<0))
%         check = 1
%         break;
%     end
% end
% E = eig(P2)
% save P2_pd P2


nav_p.xd = [0;x_2;x_3;x_4];
% nav_p.xd = [0;0;0;0];


% g = 1/((x-nav_p.xd)'*P1*(x-nav_p.xd))^(nav_p.alpha);
g = 1/((x-nav_p.xd)'*P1*(x-nav_p.xd))^(nav_p.alpha);
density = g*bump;

grad_bump = gradient(bump,[x1;x2]);

optimize = false;

matlabFunction(bump, 'File', 'functions/bump_f', 'Vars', {[x1;x2]}, 'Optimize', optimize);
matlabFunction(density, 'File', 'functions/density_f', 'Vars', {x}, 'Optimize', optimize);


matlabFunction(grad_bump, 'File', 'functions/grad_bump_f', 'Vars', {[x1;x2]}, 'Optimize', optimize);

div0 = divergence(nav_p.A*x*density,x);
matlabFunction(div0, 'File', 'functions/div0_f', 'Vars', {x}, 'Optimize', optimize);
div00 = divergence(nav_p.C*density,x);
matlabFunction(div00, 'File', 'functions/div00_f', 'Vars', {x}, 'Optimize', optimize);
div1 = divergence(nav_p.B(:,1)*density,x);
matlabFunction(div1, 'File', 'functions/div1_f', 'Vars', {x}, 'Optimize', optimize);




%% Linear programming


H_cdf = 2*eye(5*nav_p.ctrl_num+1+5);
f_cdf = zeros(5*nav_p.ctrl_num+1+5,1);
u_cdf = zeros(5*nav_p.ctrl_num+1+5,N);
x_cdf = zeros(nav_p.dim,N+1);
x_cdf_dot = zeros(nav_p.dim,N);
x_cdf(:,1) = nav_p.x0;
x1_value = zeros(nav_p.dim,1);
x2_value = zeros(nav_p.dim,1);
x3_value = zeros(nav_p.dim,1);
x4_value = zeros(nav_p.dim,1);
% x = zeros(1,N+1);
% x(1,1) = 0;
acc = zeros(1,N);
x_cdf_final = zeros(random_iter,N+1);
u_cdf_final = zeros(random_iter,N);
x_cdf_dot_final = zeros(random_iter,N);

% options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
% options = optimoptions('linprog','Algorithm','interior-point','Display','off');
options = optimoptions('quadprog','Display','off');
for rt = 1:random_iter
    
    for iter = 1:N
%         if mod(iter,10000) == 0
%             iter
%         end


        if(norm(x_cdf(:,iter)-[0;x_2;x_3;x_4])<nav_p.rad_from_goal)

            % LQR Feedback Gain
            u_cdf(1,iter) =  -K1*(x_cdf(:,iter) - nav_p.xd); 
            rf = 10;
    %         break;
    %         u_cdf(1,i) =  0;
        else
            x1_value = x_cdf(:,iter) + [epsilon;0;0;0];
            x2_value = x_cdf(:,iter) + [0;epsilon;0;0];
            x3_value = x_cdf(:,iter) + [0;0;epsilon;0];
            x4_value = x_cdf(:,iter) + [0;0;0;epsilon];


            A_cdf = [-div1_f(x_cdf(:,iter)),0,0,0,0,0,0,0,0,0,0;...
                0,-div1_f(x1_value),0,0,0,0,0,0,0,0,0;...
                0,0,-div1_f(x2_value),0,0,0,0,0,0,0,0;...
                0,0,0,-div1_f(x3_value),0,0,0,0,0,0,0;...
                0,0,0,0,-div1_f(x4_value),0,0,0,0,0,0;...
                (-sum(nav_p.B(:,1)))/epsilon,nav_p.B(1,1)/epsilon,nav_p.B(2,1)/epsilon,nav_p.B(3,1)/epsilon,nav_p.B(4,1)/epsilon,1,0,0,0,0,0;...
                (sum(nav_p.B(:,1)))/epsilon,-nav_p.B(1,1)/epsilon,-nav_p.B(2,1)/epsilon,-nav_p.B(3,1)/epsilon,-nav_p.B(4,1)/epsilon,1,0,0,0,0,0];

            b_cdf = [div0_f(x_cdf(:,iter)) - nav_p.beta*density_f(x_cdf(:,iter)) + div00_f(x_cdf(:,iter)) - kappa*density_f(x_cdf(:,iter));...
                div0_f(x1_value) - nav_p.beta*density_f(x1_value) + div00_f(x1_value) - kappa*density_f(x1_value);...
                div0_f(x2_value) - nav_p.beta*density_f(x2_value) + div00_f(x2_value) - kappa*density_f(x2_value);...
                div0_f(x3_value) - nav_p.beta*density_f(x3_value) + div00_f(x3_value) - kappa*density_f(x3_value);...
                div0_f(x4_value) - nav_p.beta*density_f(x4_value) + div00_f(x4_value) - kappa*density_f(x4_value);...
                nav_p.beta;...
                nav_p.beta];

            A_cdfe = [1,0,0,0,0,0,-1,0,0,0,0;...
                      0,1,0,0,0,0,0,-1,0,0,0;...
                      0,0,1,0,0,0,0,0,-1,0,0;...
                      0,0,0,1,0,0,0,0,0,-1,0;...
                      0,0,0,0,1,0,0,0,0,0,-1];

            B_cdfe = -[K1*(x_cdf(:,iter) - nav_p.xd);...
                       K1*(x1_value - nav_p.xd);...
                       K1*(x2_value - nav_p.xd);...
                       K1*(x3_value - nav_p.xd);...
                       K1*(x4_value - nav_p.xd)];

            max1_0 = (2*nav_p.Cf)^(-1)*(-nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x_cdf(2,iter)-(x_cdf(3,iter)*nav_p.v0)+(nav_p.a*x_cdf(4,iter)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x_cdf(2,iter)-(x_cdf(3,iter)*nav_p.v0)-(nav_p.b*x_cdf(4,iter)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x_cdf(4,iter)) - (nav_p.M*nav_p.v0*rd)) ; 
            max2_0 = (2*nav_p.Cf)^(-1)*(nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x_cdf(2,iter)-(x_cdf(3,iter)*nav_p.v0)+(nav_p.a*x_cdf(4,iter)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x_cdf(2,iter)-(x_cdf(3,iter)*nav_p.v0)-(nav_p.b*x_cdf(4,iter)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x_cdf(4,iter)) - (nav_p.M*nav_p.v0*rd));

            max1_1 = (2*nav_p.Cf)^(-1)*(-nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x1_value(2,1)-(x1_value(3,1)*nav_p.v0)+(nav_p.a*x1_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x1_value(2,1)-(x1_value(3,1)*nav_p.v0)-(nav_p.b*x1_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x1_value(4,1)) - (nav_p.M*nav_p.v0*rd)); 
            max2_1 = (2*nav_p.Cf)^(-1)*(nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x1_value(2,1)-(x1_value(3,1)*nav_p.v0)+(nav_p.a*x1_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x1_value(2,1)-(x1_value(3,1)*nav_p.v0)-(nav_p.b*x1_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x1_value(4,1)) - (nav_p.M*nav_p.v0*rd));  

            max1_2 = (2*nav_p.Cf)^(-1)*(-nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x2_value(2,1)-(x2_value(3,1)*nav_p.v0)+(nav_p.a*x2_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x2_value(2,1)-(x2_value(3,1)*nav_p.v0)-(nav_p.b*x2_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x2_value(4,1)) - (nav_p.M*nav_p.v0*rd)); 
            max2_2 = (2*nav_p.Cf)^(-1)*(nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x2_value(2,1)-(x2_value(3,1)*nav_p.v0)+(nav_p.a*x2_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x2_value(2,1)-(x2_value(3,1)*nav_p.v0)-(nav_p.b*x2_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x2_value(4,1)) - (nav_p.M*nav_p.v0*rd));

            max1_3 = (2*nav_p.Cf)^(-1)*(-nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x3_value(2,1)-(x3_value(3,1)*nav_p.v0)+(nav_p.a*x3_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x3_value(2,1)-(x3_value(3,1)*nav_p.v0)-(nav_p.b*x3_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x3_value(4,1)) - (nav_p.M*nav_p.v0*rd)); 
            max2_3 = (2*nav_p.Cf)^(-1)*(nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x3_value(2,1)-(x3_value(3,1)*nav_p.v0)+(nav_p.a*x3_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x3_value(2,1)-(x3_value(3,1)*nav_p.v0)-(nav_p.b*x3_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x3_value(4,1)) - (nav_p.M*nav_p.v0*rd));

            max1_4 = (2*nav_p.Cf)^(-1)*(-nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x4_value(2,1)-(x4_value(3,1)*nav_p.v0)+(nav_p.a*x4_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x4_value(2,1)-(x4_value(3,1)*nav_p.v0)-(nav_p.b*x4_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x4_value(4,1)) - (nav_p.M*nav_p.v0*rd)); 
            max2_4 = (2*nav_p.Cf)^(-1)*(nav_p.M*nav_p.acc_max + (2*nav_p.Cf*(x4_value(2,1)-(x4_value(3,1)*nav_p.v0)+(nav_p.a*x4_value(4,1)))/(nav_p.v0)) ...
                                    + (2*(nav_p.Cr+Cr_max)*(x4_value(2,1)-(x4_value(3,1)*nav_p.v0)-(nav_p.b*x4_value(4,1)))/(nav_p.v0)) + (2*nav_p.M*nav_p.v0*x4_value(4,1)) - (nav_p.M*nav_p.v0*rd));





            lb_cdf = [max1_0;max1_1;max1_2;max1_3;max1_4;1e-10;-delta_bound;-delta_bound;-delta_bound;-delta_bound;-delta_bound];
            ub_cdf = [max2_0;max2_1;max2_2;max2_3;max2_4;Inf;delta_bound;delta_bound;delta_bound;delta_bound;delta_bound];


            u_cdf(:,iter) = ctrl_multiplier*quadprog(H_cdf,f_cdf,A_cdf,b_cdf,[],[],lb_cdf,ub_cdf,[],options);


        end
        nav_p.A_uncertain = [0,1,0,-nav_p.L;...
           0,-(2*((nav_p.Cr+rand_Cr(rt,1))+nav_p.Cf))/(nav_p.M*nav_p.v0),(2*((nav_p.Cr+rand_Cr(rt,1))+nav_p.Cf))/(nav_p.M),((2*((nav_p.b*(nav_p.Cr+rand_Cr(rt,1))) - (nav_p.a*nav_p.Cf)))/(nav_p.M*nav_p.v0))-(2*nav_p.v0);...
           0,0,0,-1;...
           0,(2*((nav_p.b*(nav_p.Cr+rand_Cr(rt,1))) - (nav_p.a*nav_p.Cf)))/(nav_p.Iz*nav_p.v0),(-2*((nav_p.b*(nav_p.Cr+rand_Cr(rt,1))) - (nav_p.a*nav_p.Cf)))/(nav_p.Iz),(-2*((nav_p.a^2)*nav_p.Cf + (nav_p.b^2)*(nav_p.Cr+rand_Cr(rt,1))))/(nav_p.Iz*nav_p.v0)];
        
        x_cdf_dot(:,iter) = nav_p.A_uncertain*x_cdf(:,iter) + nav_p.B*u_cdf(1,iter) + nav_p.C;
        x_cdf(:,iter+1) = x_cdf(:,iter) + (nav_p.A_uncertain*x_cdf(:,iter) + nav_p.B*u_cdf(1,iter)+ nav_p.C)*deltaT;


    end
    x_cdf_final(rt,:) = x_cdf(1,:);
    u_cdf_final(rt,:) = u_cdf(1,:);
    x_cdf_dot_final(rt,:) = x_cdf_dot(2,:);    
end
%% Plot


% t = 1:size(u_cdf,2);
% t1 = 1:size(x_cdf,2);
t = 1:iter;
t1 = 1:iter+1;


figure()
plot(t1*deltaT,x_cdf_final(:,1:N+1),'Color',[0.2 0.5 0.9 alpha], 'LineWidth', 2); hold on;
plot(t1*deltaT,nav_p.disp_max_h*ones(1,size(t1,2)),'r', 'LineWidth', 2); hold on;
plot(t1*deltaT,nav_p.disp_max_s*ones(1,size(t1,2)),'g', 'LineWidth', 2); hold on;
plot(t1*deltaT,-nav_p.disp_max_h*ones(1,size(t1,2)),'r', 'LineWidth', 2); hold on;
plot(t1*deltaT,-nav_p.disp_max_s*ones(1,size(t1,2)),'g', 'LineWidth', 2); hold on;
hold off;
title("Lateral Displacement")
xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$x_1$",'interpreter','latex', 'FontSize', 20)
ctrl_lgd = legend('$x_1$','r1 bound','r2 bound','interpreter','latex');
ctrl_lgd.FontSize = 20;
xlim([0,N*deltaT]);
ylim([-nav_p.disp_max_h-0.1,nav_p.disp_max_h+0.1]);
grid on;

figure()
plot(t*deltaT,u_cdf_final(:,1:N),'Color',[0.2 0.5 0.9 alpha], 'LineWidth', 2); hold on;
hold off;
title("control plot")
xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$u$",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend('acc [g]','bounds','interpreter','latex');
% ctrl_lgd.FontSize = 20;
xlim([0,N*deltaT]);
grid on;

figure()
plot(t*deltaT,x_cdf_dot_final(:,1:N)/9.81,'Color',[0.2 0.5 0.9 alpha], 'LineWidth', 2); hold on;
plot(t*deltaT,nav_p.acc_max*ones(1,size(t,2))/9.81,'r', 'LineWidth', 2); hold on;
plot(t*deltaT,-nav_p.acc_max*ones(1,size(t,2))/9.81,'r', 'LineWidth', 2); hold on;
hold off;
title("Lateral Accelaration [g]")
xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$\ddot{x_1}$ [g]",'interpreter','latex', 'FontSize', 20)
ctrl_lgd = legend('acc [g]','bounds','interpreter','latex');
ctrl_lgd.FontSize = 20;
xlim([0,N*deltaT]);
ylim([-(nav_p.acc_max/9.81)-0.1,(nav_p.acc_max/9.81)+0.1]);
grid on;

figure()
plot(t1*deltaT,x_cdf(1,1:iter+1), 'LineWidth', 2); hold on;
plot(t1*deltaT,x_cdf(2,1:iter+1), 'LineWidth', 2); hold on;
plot(t1*deltaT,x_cdf(3,1:iter+1), 'LineWidth', 2); hold on;
plot(t1*deltaT,x_cdf(4,1:iter+1), 'LineWidth', 2); hold on;
hold off;
title("state trajectory")
xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$states$",'interpreter','latex', 'FontSize', 20)
ctrl_lgd = legend('$x_1$','$x_2$','$x_3$','$x_4$','interpreter','latex');
ctrl_lgd.FontSize = 14;
xlim([0,N*deltaT]);
% ylim([-20,20]);
grid on;

%% Plotting Density & Rantzer Conditions
x_vec = -25:0.1:25; y_vec = x_vec;
[X,Y] = meshgrid(x_vec,y_vec);
rho_val = zeros(size(X));
grad_bump_x1 = zeros(size(X));
grad_bump_x2 = zeros(size(Y));

for i=1:length(x_vec)
    for j = 1:length(y_vec)
        rho_val(i,j) = bump_f([x_vec(j);y_vec(i)]);
        grad_bump = grad_bump_f([x_vec(j);y_vec(i)]);
        grad_bump_x1(i,j) = grad_bump(1);
        grad_bump_x2(i,j) = grad_bump(2);
    end
end
figure()
surf(X,Y,rho_val, 'FaceAlpha',0.65, 'EdgeColor', 'none')
colormap jet
view(90,60)
title("Density Function")
xlabel("x_1")
ylabel("x_2")
zlim([-1 3])

figure()
quiverInLogScale(X, Y, grad_bump_x1, grad_bump_x2); % Add 3 for scaling
title("Log Scale Gradient of Density Function");
xlabel("x_1")
ylabel("x_2")
hold on;
