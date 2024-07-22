function bump = formOval(r1, r2, a, b, c, k, x, p, vpa_enable, theta)
%formCircleBump
% Forms circular bump function of R^n based off obstacle radius (r1), sensing
% radius (epsilon = r2-r1), center, and MATLAB symbolic variable x
%
% Inputs:
%   r1          : Positive value of radius of obstacle (R^1)
%   r2          : Positive value of obstacle radius + sensing radius
%                   r2 = r1 + epsilon (R^1)
%   c           : Center of bump function (R^n)
%   x           : Symbolic variable of R^n (n determines dimension of bump)
%   p           : Scalar to choose the p-norm
%   vpa_enable  : Boolean to enable VPA precision w/ digit 10
%   A           : Diffeomorphic mapping (R^(nxn))
%
% Outputs:
%   bump        : Symbolic expression of bump

% TODO: Figure out how to use compose s.t. mathematically cleaner...
% currently having trouble w/ replicating results from compose
if nargin < 5
    p = 2; % Default 2-norm
end
if nargin < 6
    vpa_enable = false;
end
if (nargin < 7)
    theta = 0;
end




% Linear change of coordinates and symmetric w/ norm(x)^2

term1 = ((x(1)-c(1))*cos(theta)-(x(2)-c(2))*sin(theta))^2/(a^2);
term2 = ((x(1)-c(1))*sin(theta)+(x(2)-c(2))*cos(theta))^2/(b^2);
% term3 = (1+k*((x(1)-c(1))*cos(theta)-(x(2)-c(2))*sin(theta)))/(1-k*((x(1)-c(1))*cos(theta)-(x(2)-c(2))*sin(theta)));
term3 = k^((x(1)-c(1))*cos(theta)-(x(2)-c(2))*sin(theta));

mm = term1 + (term2*term3)-r1;

m = mm/(r2-r1);

f = piecewise(m<=0, 0, exp(-1/m));
f_shift = piecewise(1-m <= 0, 0, exp(-1/(1-m)));

% Normalize the symmetric step function
g = f/(f+f_shift);
bump = g; % Emphasis that g is bump function

if vpa_enable
    bump = vpa(bump,5);
end


end