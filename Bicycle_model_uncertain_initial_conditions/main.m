%% Obstacle Avoidance control using quadratic program
clc; clear; close all; 
mkdir('functions');
addpath('./functions');

%% Initialization
nav_p.x0 = [-6;0];%[-6;0.1;0;0]; % Start position
nav_p.xd = [6;-0.5];%[6;0;0;0]; % Goal position / Desired position
nav_p.p = 2; % Generate obstacle set off p-norm
nav_p.rad_from_goal = 0.01;
num_random_init = 20;
max_rand = 0.5;
rng(2);
% random_init = nav_p.x0 + (2*max_rand*rand(2,num_random_init)) - max_rand;
random_init = [-6.064,-5.9503,-6.0796,-6.2954,-5.9503,-5.9503,-6.0796,-6.0796,-6.064,-6.064,-6.2954,-6.2954,-6.2954,-6.2954,-6.0796,-6.0796,-6.0796,-6.2954,-6.2954,-6.2954;...
               -0.4741,-0.0647,-0.1697,0.1193,-0.0947,-0.1347,-0.2097,-0.2597,-0.4241,-0.3841,0.0793,0.0393,0.1693,0.2093,-0.2997,-0.3297,-0.3697,0.2493,0.2893,0.3293];
nav_p.r1 = 1.3;
nav_p.r2 = 1.8;
nav_p.alpha = 0.2;
nav_p.ctrl_num = 2;
nav_p.dim = 2;
nav_p.A = zeros(2);
nav_p.B = eye(2);
nav_p.Q = eye(2);
nav_p.R = 1;
nav_p.alpha1 = 2;
nav_p.alpha2 = 30;

lr = 1.2;
L = 2;

N = 0.8e3; %timesteps
deltaT = 1e-2;%1e-2;
ctrl_bound = 1e5;
a_ctrl_bound = 1e5;
w_ctrl_bound = 1e5;
theta = 0*pi/180; % [Radian] CounterClockwise | Rotate after
stretch = [1;1];
epsilon = 0.001;


% Obstacle centers
num_obs = 2;

nav_p.c1 = [-2; 0.5]; 
nav_p.c2 = [1.7; -0.8];


%% Density Function Formulation
syms x [nav_p.dim,1] real

[A1, A_inv] = transformationMatrix(theta, stretch, 2);

[K,P1,e] = lqr(nav_p.A,nav_p.B,nav_p.Q,nav_p.R);


bump = formPNormBump(nav_p.r1,nav_p.r2, nav_p.c1, x, nav_p.p, true)*...
    formPNormBump(nav_p.r1,nav_p.r2, nav_p.c2, [x1;x2], nav_p.p, true);%*...
%     formPNormBump(nav_p.r1,nav_p.r2, nav_p.c3, [x1;x2], nav_p.p, true);

% bump = formFastInvBump([x1;x2], nav_p, num_obs);


g = 1/((x-nav_p.xd)'*P1*(x-nav_p.xd))^(nav_p.alpha);
density = g*bump;
grad_density = gradient(density, x);
hess_density = hessian(density, x);

optimize = false;

matlabFunction(bump, 'File', 'functions/bump_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(density, 'File', 'functions/density_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(grad_density, 'File', 'functions/grad_density_f', 'Vars', {x}, 'Optimize', optimize);
matlabFunction(hess_density, 'File', 'functions/hess_density_f', 'Vars', {x}, 'Optimize', optimize);


div0 = divergence(nav_p.A*x*density,x);
matlabFunction(div0, 'File', 'functions/div0_f', 'Vars', {x}, 'Optimize', optimize);
div1 = divergence(nav_p.B(:,1)*density,x);
matlabFunction(div1, 'File', 'functions/div1_f', 'Vars', {x}, 'Optimize', optimize);
div2 = divergence(nav_p.B(:,2)*density,x);
matlabFunction(div2, 'File', 'functions/div2_f', 'Vars', {x}, 'Optimize', optimize);



%% Linear programming
decision_varaible = (nav_p.dim+1)*nav_p.ctrl_num+1; % Right now no 'c' slack
                                                    % 's' is slack to make rantzer condition positive
                                                    % 'c' is required to make constraint work if density = bump 
H_cdf = 2*blkdiag(1,1,0,0,0,0,1);
% H_cdf = 2*blkdiag(1,1,1,1,1,1,1);
f_cdf = zeros(decision_varaible,1);
u_cdf = NaN*ones(decision_varaible,N);
u_cdf_certain = NaN*ones(decision_varaible,N);
x_cdf = NaN*ones(nav_p.dim,N+1);
x_cdf_certain = NaN*ones(nav_p.dim,N+1);
x_cdf(:,1) = nav_p.x0;
x_cdf_certain(:,1) = nav_p.x0;

x_cdf_extra = zeros(nav_p.dim*num_random_init,N+1);
x_cdf_extra(:,1) = random_init(:);

theta = NaN*ones(1,N+1);
theta_certain = NaN*ones(1,N+1);
theta(:,1) = 0*(pi/180);
theta_certain(:,1) = 0*(pi/180);

velocity_lin = NaN*ones(1,N+1);
velocity_lin_certain = NaN*ones(1,N+1);
velocity_lin(:,1) = 0;
velocity_lin_certain(:,1) = 0;

steering_angle = NaN*ones(1,N+1);
steering_angle_certain = NaN*ones(1,N+1);
steering_angle(:,1) = 0;
steering_angle_certain(:,1) = 0;

v_tilda = NaN*ones(1,N);
v_tilda_certain = NaN*ones(1,N);

v_tilda_dot = NaN*ones(1,N);
v_tilda_dot_certain = NaN*ones(1,N);

Phi_tilda = NaN*ones(1,N);
Phi_tilda_certain = NaN*ones(1,N);
Phi_tilda_dot = NaN*ones(1,N);
Phi_tilda_dot_certain = NaN*ones(1,N);

w_x = NaN*ones(1,N);
w_x_certain = NaN*ones(1,N);
a_x = NaN*ones(1,N);
a_x_certain = NaN*ones(1,N);


x1_value = NaN*ones(nav_p.dim,1);
x2_value = NaN*ones(nav_p.dim,1);

x1_value_certain = NaN*ones(nav_p.dim,1);
x2_value_certain = NaN*ones(nav_p.dim,1);

x1_value_extra = zeros(nav_p.dim*num_random_init,1);
x2_value_extra = zeros(nav_p.dim*num_random_init,1);


% options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
options = optimoptions('quadprog','Display','off','MaxIterations',1e6,'Algorithm','interior-point-convex');

for i = 1:N
    if mod(i,1e4) == 0
        i
    end
    
    if(norm(x_cdf(:,i)-nav_p.xd)<nav_p.rad_from_goal)
 
        w_x(:,i) = 0;
        a_x(:,i) = 0;
        velocity_lin(:,i) = 0;
    
    else
    
      
        x1_value = x_cdf(:,i) + [epsilon;0];
        x2_value = x_cdf(:,i) + [0;epsilon]; 

        x1_value_extra = x_cdf_extra(:,i) + repmat([epsilon;0],num_random_init,1);
        x2_value_extra = x_cdf_extra(:,i) + repmat([0;epsilon],num_random_init,1);        

        A_cdf = [-div1_f(x_cdf(:,i)),-div2_f(x_cdf(:,i)),0,0,0,0,1;...
                 0,0,-div1_f(x1_value),-div2_f(x1_value),0,0,1;...
                 0,0,0,0,-div1_f(x2_value),-div2_f(x2_value),1;...
                 sum(nav_p.B(:,1))/epsilon,sum(nav_p.B(:,2))/epsilon,-nav_p.B(1,1)/epsilon,-nav_p.B(1,2)/epsilon,-nav_p.B(2,1)/epsilon,-nav_p.B(2,2)/epsilon,0;...
                 DIV_A_original(x_cdf_extra(:,i),nav_p.dim,num_random_init);...
                 DIV_A_x1(x1_value_extra,nav_p.dim,num_random_init);...
                 DIV_A_x2(x2_value_extra,nav_p.dim,num_random_init)];

        b_cdf = [div0_f(x_cdf(:,i));div0_f(x1_value);div0_f(x2_value);0;...
                 DIV_B(x_cdf_extra(:,i),x1_value_extra,x2_value_extra,nav_p.dim,num_random_init,nav_p.A)];

        lb_cdf = [-ctrl_bound*ones(6,1);1e-5];
        ub_cdf = [ctrl_bound*ones(6,1);1e5];

        feedback_x = -K*(x_cdf(:,i)-nav_p.xd);

        f_cdf = -2*[feedback_x;zeros(4,1);0];

        u_cdf(:,i) = quadprog(H_cdf,f_cdf,A_cdf,b_cdf,[],[],lb_cdf,ub_cdf,[],options);
        

        v_tilda(1,i) = sqrt(u_cdf(1,i)^2+u_cdf(2,i)^2);

        if i == 1
            v_tilda_dot(1,i) = (v_tilda(1,i)-0)/deltaT;
        else
            v_tilda_dot(1,i) = (v_tilda(1,i)-v_tilda(1,i-1))/deltaT;
        end 

        a_x(1,i) = v_tilda_dot(1,i) - nav_p.alpha1*(velocity_lin(1,i)-v_tilda(1,i));
        
        
        if abs(a_x(1,i)) >= a_ctrl_bound
            a_x(1,i) = sign(a_x(1,i))*a_ctrl_bound;
        end

        Phi_tilda(1,i) = atan2(u_cdf(2,i),u_cdf(1,i)); % value wrapped between -pi and pi


        if i == 1
            Phi_tilda_dot(1,i) = (Phi_tilda(1,i)-0)/deltaT;
        else
            Phi_tilda_dot(1,i) = (Phi_tilda(1,i)-Phi_tilda(1,i-1))/deltaT;
        end    

        w_x(1,i) = ((1+((lr*tan(steering_angle(:,i)))/(L))^2)/((lr*(sec(steering_angle(:,i)))^2)/(L)))*...
                   ((-velocity_lin(:,i))/(L)*(cos(atan2(lr*tan(steering_angle(:,i)),L))*tan(steering_angle(:,i))) + Phi_tilda_dot(1,i) - nav_p.alpha2*sin(theta(:,i) + atan2(lr*tan(steering_angle(:,i)),L) - Phi_tilda(1,i)));

       
        if abs(w_x(1,i)) >= w_ctrl_bound
            w_x(1,i) = sign(w_x(1,i))*w_ctrl_bound;
        end
        
    end
    x_cdf(:,i+1) = x_cdf(:,i) + ([velocity_lin(:,i)*cos(theta(:,i)+atan2(lr*tan(steering_angle(:,i)),L));velocity_lin(:,i)*sin(theta(:,i)+atan2(lr*tan(steering_angle(:,i)),L))])*deltaT;
    x_cdf_extra(:,i+1) = x_cdf_extra(:,i) + repmat(([velocity_lin(:,i)*cos(theta(:,i)+atan2(lr*tan(steering_angle(:,i)),L));velocity_lin(:,i)*sin(theta(:,i)+atan2(lr*tan(steering_angle(:,i)),L))])*deltaT,num_random_init,1);
    theta(:,i+1) = theta(:,i) + (velocity_lin(:,i)*cos(atan2(lr*tan(steering_angle(:,i)),L))*tan(steering_angle(:,i))/L)*deltaT;
    steering_angle(:,i+1) = steering_angle(:,i) + (w_x(1,i))*deltaT;
    velocity_lin(:,i+1) = velocity_lin(:,i) + (a_x(1,i))*deltaT;
    
%     theta(:,i+1) = pi - mod(theta(:,i+1),2*pi); %wrapping between -pi and pi
    if theta(1,i+1) > pi
        j = 1;
        while theta(1,i+1)/pi > (2*j) + 1
            j = j+1;
        end
        theta(1,i+1) = theta(1,i+1) - j*(2*pi);
    end
    if theta(1,i+1) < -pi
        j = 1;
        while theta(1,i+1)/(-pi) > (2*j) + 1
            j = j+1;
        end
        theta(1,i+1) = theta(1,i+1) + j*(2*pi);
    end
    
    if(norm(x_cdf_certain(:,i)-nav_p.xd)<nav_p.rad_from_goal)
 
        w_x_certain(:,i) = 0;
        a_x_certain(:,i) = 0;
        velocity_lin_certain(:,i) = 0;
    
    else
    
      
        x1_value_certain = x_cdf_certain(:,i) + [epsilon;0];
        x2_value_certain = x_cdf_certain(:,i) + [0;epsilon]; 
        

        A_cdf_certain = [-div1_f(x_cdf_certain(:,i)),-div2_f(x_cdf_certain(:,i)),0,0,0,0,1;...
                 0,0,-div1_f(x1_value_certain),-div2_f(x1_value_certain),0,0,1;...
                 0,0,0,0,-div1_f(x2_value_certain),-div2_f(x2_value_certain),1;...
                 sum(nav_p.B(:,1))/epsilon,sum(nav_p.B(:,2))/epsilon,-nav_p.B(1,1)/epsilon,-nav_p.B(1,2)/epsilon,-nav_p.B(2,1)/epsilon,-nav_p.B(2,2)/epsilon,0];

        b_cdf_certain = [div0_f(x_cdf_certain(:,i));div0_f(x1_value_certain);div0_f(x2_value_certain);0];

        lb_cdf = [-ctrl_bound*ones(6,1);1e-5];
        ub_cdf = [ctrl_bound*ones(6,1);1e5];

        feedback_x_certain = -K*(x_cdf_certain(:,i)-nav_p.xd);

        f_cdf = -2*[feedback_x_certain;zeros(4,1);0];

        u_cdf_certain(:,i) = quadprog(H_cdf,f_cdf,A_cdf_certain,b_cdf_certain,[],[],lb_cdf,ub_cdf,[],options);
        

        v_tilda_certain(1,i) = sqrt(u_cdf_certain(1,i)^2+u_cdf_certain(2,i)^2);

        if i == 1
            v_tilda_dot_certain(1,i) = (v_tilda_certain(1,i)-0)/deltaT;
        else
            v_tilda_dot_certain(1,i) = (v_tilda_certain(1,i)-v_tilda_certain(1,i-1))/deltaT;
        end 

        a_x_certain(1,i) = v_tilda_dot_certain(1,i) - nav_p.alpha1*(velocity_lin_certain(1,i)-v_tilda_certain(1,i));
        
        
        if abs(a_x_certain(1,i)) >= a_ctrl_bound
            a_x_certain(1,i) = sign(a_x_certain(1,i))*a_ctrl_bound;
        end

        Phi_tilda_certain(1,i) = atan2(u_cdf_certain(2,i),u_cdf_certain(1,i)); % value wrapped between -pi and pi


        if i == 1
            Phi_tilda_dot_certain(1,i) = (Phi_tilda_certain(1,i)-0)/deltaT;
        else
            Phi_tilda_dot_certain(1,i) = (Phi_tilda_certain(1,i)-Phi_tilda_certain(1,i-1))/deltaT;
        end    

        w_x_certain(1,i) = ((1+((lr*tan(steering_angle_certain(:,i)))/(L))^2)/((lr*(sec(steering_angle_certain(:,i)))^2)/(L)))*...
                   ((-velocity_lin_certain(:,i))/(L)*(cos(atan2(lr*tan(steering_angle_certain(:,i)),L))*tan(steering_angle_certain(:,i))) + Phi_tilda_dot_certain(1,i) - nav_p.alpha2*sin(theta_certain(:,i) + atan2(lr*tan(steering_angle_certain(:,i)),L) - Phi_tilda_certain(1,i)));

       
        if abs(w_x_certain(1,i)) >= w_ctrl_bound
            w_x_certain(1,i) = sign(w_x_certain(1,i))*w_ctrl_bound;
        end
        
    end
    x_cdf_certain(:,i+1) = x_cdf_certain(:,i) + ([velocity_lin_certain(:,i)*cos(theta_certain(:,i)+atan2(lr*tan(steering_angle_certain(:,i)),L));velocity_lin_certain(:,i)*sin(theta_certain(:,i)+atan2(lr*tan(steering_angle_certain(:,i)),L))])*deltaT;
    theta_certain(:,i+1) = theta_certain(:,i) + (velocity_lin_certain(:,i)*cos(atan2(lr*tan(steering_angle_certain(:,i)),L))*tan(steering_angle_certain(:,i))/L)*deltaT;
    steering_angle_certain(:,i+1) = steering_angle_certain(:,i) + (w_x_certain(1,i))*deltaT;
    velocity_lin_certain(:,i+1) = velocity_lin_certain(:,i) + (a_x_certain(1,i))*deltaT;
    
%     theta(:,i+1) = pi - mod(theta(:,i+1),2*pi); %wrapping between -pi and pi
    if theta_certain(1,i+1) > pi
        j = 1;
        while theta_certain(1,i+1)/pi > (2*j) + 1
            j = j+1;
        end
        theta_certain(1,i+1) = theta_certain(1,i+1) - j*(2*pi);
    end
    if theta_certain(1,i+1) < -pi
        j = 1;
        while theta_certain(1,i+1)/(-pi) > (2*j) + 1
            j = j+1;
        end
        theta_certain(1,i+1) = theta_certain(1,i+1) + j*(2*pi);
    end
    

end

%% Plot
close all;

obsColor = [.7 .7 .7]; % Obstacle color -> Grey
alpha = 0.2;

x_vec = -11:0.25:11; y_vec = x_vec;
[X,Y] = meshgrid(x_vec,y_vec);
sensing_radius_c1 = zeros(size(X));
sensing_radius_c2 = zeros(size(X));
uncertain_radius = zeros(size(X));
target_radius = zeros(size(X));

for i=1:length(x_vec)
    for j = 1:length(y_vec)
        sensing_radius_c1(i,j) = norm(eye(2)*([x_vec(j);y_vec(i)]-nav_p.c1(1:2,:)), nav_p.p)-nav_p.r2;
        sensing_radius_c2(i,j) = norm(eye(2)*([x_vec(j);y_vec(i)]-nav_p.c2(1:2,:)), nav_p.p)-nav_p.r2;
        uncertain_radius(i,j) = norm(eye(2)*([x_vec(j);y_vec(i)]-nav_p.x0), nav_p.p)-max_rand;
        target_radius(i,j) = norm(eye(2)*([x_vec(j);y_vec(i)]-nav_p.xd), nav_p.p)-max_rand;

    end
end

figure()
plot(x_cdf(1,:),x_cdf(2,:),'Color',[0.2 0.5 0.9 alpha], 'LineWidth', 2); hold on;
plot(nav_p.x0(1), nav_p.x0(2), 'ok', 'MarkerSize',10, 'MarkerFaceColor',"#D95319"); hold on;
plot(nav_p.xd(1), nav_p.xd(2), 'ok', 'MarkerSize',10, 'MarkerFaceColor',"#77AC30"); hold on;
plot(NaN,NaN, 'o','MarkerSize', 10, 'MarkerEdgeColor',...
            'black', 'MarkerFaceColor',obsColor, 'LineWidth', 1.5); % For legend as rectangular object can't be defined as a legend
      
contour(X,Y,sensing_radius_c1,[0 0],'--k','LineWidth',1); hold on;  

gamma = (0:100-1)*(2*pi/100);
points = [nav_p.c1(1);nav_p.c1(2)] + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;
points = nav_p.c2 + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;


contour(X,Y,sensing_radius_c2,[0 0],'--k','LineWidth',1); hold on;
contour(X,Y,uncertain_radius,[0 0],'Color','#D95319','LineStyle','--','LineWidth',1); hold on;
contour(X,Y,target_radius,[0 0],'Color',"#77AC30",'LineStyle','--','LineWidth',1); hold on 

total_traj_success = 0;
for i = 1:num_random_init
    if(norm(x_cdf_extra((nav_p.dim*i) - (nav_p.dim-1):(nav_p.dim*i) - (nav_p.dim-1) + 1,N)-nav_p.xd)<max_rand)
        total_traj_success = total_traj_success + 1;
        plot(x_cdf_extra((nav_p.dim*i) - (nav_p.dim-1),:),x_cdf_extra((nav_p.dim*i) - (nav_p.dim-1) + 1,:),'Color',[0.2 0.5 0.9 alpha], 'LineWidth', 2); hold on;
    end
end
grid off;
hold off;
title("Safe navigation with uncertainty",'interpreter','latex', 'FontSize', 35)
xlabel('position $x_1$ [m]','interpreter','latex', 'FontSize', 30);
ylabel('position $x_2$ [m]','interpreter','latex', 'FontSize', 30);
xlim([-6,6]);
ylim([-6,6]);
axis square

figure()
plot(x_cdf_certain(1,:),x_cdf_certain(2,:),'Color','#D95319', 'LineStyle','--','LineWidth', 2); hold on;
plot(nav_p.x0(1), nav_p.x0(2), 'ok', 'MarkerSize',10, 'MarkerFaceColor',"#D95319"); hold on;
plot(nav_p.xd(1), nav_p.xd(2), 'ok', 'MarkerSize',10, 'MarkerFaceColor',"#77AC30"); hold on;
plot(NaN,NaN, 'o','MarkerSize', 10, 'MarkerEdgeColor',...
            'black', 'MarkerFaceColor',obsColor, 'LineWidth', 1.5); % For legend as rectangular object can't be defined as a legend
      
contour(X,Y,sensing_radius_c1,[0 0],'--k','LineWidth',1); hold on;  

gamma = (0:100-1)*(2*pi/100);
points = [nav_p.c1(1);nav_p.c1(2)] + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;
points = nav_p.c2 + A1*[nav_p.r1*cos(gamma);nav_p.r1*sin(gamma)];
P = polyshape(points(1,:), points(2,:));
plot(P, 'FaceColor', obsColor, 'LineWidth', 2, 'FaceAlpha', 1.0); hold on;


contour(X,Y,sensing_radius_c2,[0 0],'--k','LineWidth',1); hold on;


grid off;
hold off;
title("Safe navigation without uncertainty",'interpreter','latex', 'FontSize', 35)
xlabel('position $x_1$ [m]','interpreter','latex', 'FontSize', 30);
ylabel('position $x_2$ [m]','interpreter','latex', 'FontSize', 30);
xlim([-6,6]);
ylim([-6,6]);
axis square

t = 1:size(u_cdf,2);
t1 = 1:size(x_cdf,2);

figure()
subplot(3,1,1)
plot(t1*deltaT,theta(1,:),'Color','#0072BD','LineStyle','-', 'LineWidth', 2); hold on;
plot(t1*deltaT,theta_certain(1,:),'Color','#D95319', 'LineWidth', 2); hold on;
hold off;
% xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$\theta$ [rad]",'interpreter','latex', 'FontSize', 45)
ctrl_lgd = legend('uncertain','certain','interpreter','latex');
ctrl_lgd.FontSize = 20;
xlim([0,N*deltaT]);
grid off;

subplot(3,1,2)
plot(t1*deltaT,steering_angle(1,:),'Color','#0072BD','LineStyle','-', 'LineWidth', 2); hold on;
plot(t1*deltaT,steering_angle_certain(1,:),'Color','#D95319', 'LineWidth', 2); hold on;
hold off;
% xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
ylabel("$\delta$ [rad]",'interpreter','latex', 'FontSize', 45)
ctrl_lgd = legend('uncertain','certain','interpreter','latex');
ctrl_lgd.FontSize = 20;
xlim([0,N*deltaT]);
grid off;

subplot(3,1,3)
plot(t1*deltaT,velocity_lin(1,:),'Color','#0072BD','LineStyle','-', 'LineWidth', 2); hold on;
plot(t1*deltaT,velocity_lin_certain(1,:),'Color','#D95319', 'LineWidth', 2); hold on;
hold off;
xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 45)
ylabel("$v$ [mps]",'interpreter','latex', 'FontSize', 45)
ctrl_lgd = legend('uncertain','certain','interpreter','latex');
ctrl_lgd.FontSize = 20;
xlim([0,N*deltaT]);
grid off;

% 
% figure()
% plot(t*deltaT,w_x(1,:),'b', 'LineWidth', 2); hold on;
% plot(t*deltaT,a_x(1,:),'r', 'LineWidth', 2); hold on;
% hold off;
% % title("u1 Control trajectory")
% xlabel("$time$ [s]",'interpreter','latex', 'FontSize', 20)
% ylabel("control",'interpreter','latex', 'FontSize', 20)
% ctrl_lgd = legend("$\omega$","$acc$",'interpreter','latex');
% ctrl_lgd.FontSize = 14;
% xlim([0,N*deltaT]);
% grid on;

%% Block Functions

function yr = DIV_A_original(xr,dim,num_random_init)

DIV1_F = zeros(num_random_init,1);
DIV2_F = zeros(num_random_init,1);

for i = 1:num_random_init
    DIV1_F(i,1) = div1_f(xr((dim*i) - (dim-1):(dim*i),1));
    DIV2_F(i,1) = div2_f(xr((dim*i) - (dim-1):(dim*i),1));
    
end

yr = [-DIV1_F,-DIV2_F,zeros(num_random_init,4),ones(num_random_init,1)];

end

function yr = DIV_A_x1(xr,dim,num_random_init)

DIV1_F = zeros(num_random_init,1);
DIV2_F = zeros(num_random_init,1);

for i = 1:num_random_init
    DIV1_F(i,1) = div1_f(xr((dim*i) - (dim-1):(dim*i),1));
    DIV2_F(i,1) = div2_f(xr((dim*i) - (dim-1):(dim*i),1));
    
end

yr = [zeros(num_random_init,2),-DIV1_F,-DIV2_F,zeros(num_random_init,2),ones(num_random_init,1)];

end

function yr = DIV_A_x2(xr,dim,num_random_init)

DIV1_F = zeros(num_random_init,1);
DIV2_F = zeros(num_random_init,1);

for i = 1:num_random_init
    DIV1_F(i,1) = div1_f(xr((dim*i) - (dim-1):(dim*i),1));
    DIV2_F(i,1) = div2_f(xr((dim*i) - (dim-1):(dim*i),1));
    
end

yr = [zeros(num_random_init,4),-DIV1_F,-DIV2_F,ones(num_random_init,1)];

end

function yr = DIV_B(xr,xr1,xr2,dim,num_random_init,A)

p = zeros(num_random_init,1);
p1 = zeros(num_random_init,1);
p2 = zeros(num_random_init,1);

for i = 1:num_random_init
    p(i,1) = trace(A)*density_f(xr((dim*i) - (dim-1):(dim*i),1)) + grad_density_f(xr((dim*i) - (dim-1):(dim*i),1))'*A*xr((dim*i) - (dim-1):(dim*i),1);
    p1(i,1) = trace(A)*density_f(xr1((dim*i) - (dim-1):(dim*i),1)) + grad_density_f(xr1((dim*i) - (dim-1):(dim*i),1))'*A*xr1((dim*i) - (dim-1):(dim*i),1);
    p2(i,1) = trace(A)*density_f(xr2((dim*i) - (dim-1):(dim*i),1)) + grad_density_f(xr2((dim*i) - (dim-1):(dim*i),1))'*A*xr2((dim*i) - (dim-1):(dim*i),1);

    
end

yr = [p;p1;p2];

end