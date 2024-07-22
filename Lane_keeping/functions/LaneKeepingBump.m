function bump = LaneKeepingBump(r1, r2, a, x)
%formCircleBump
% Forms circular bump function of R^n based off obstacle radius (r1), sensing
% radius (epsilon = r2-r1), center, and MATLAB symbolic variable x
%
% Inputs:
%   r1          : Positive value of radius of obstacle (R^1)
%   r2          : Positive value of obstacle radius + sensing radius
%                   r2 = r1 + epsilon (R^1)
%   x           : Symbolic variable of R^n (n determines dimension of bump)
%   a           : Maximum acceleration value

%
% Outputs:
%   bump        : Symbolic expression of bump

% TODO: Figure out how to use compose s.t. mathematically cleaner...
% currently having trouble w/ replicating results from compose

if nargin < 6
    vpa_enable = false;
end





% Linear change of coordinates and symmetric w/ norm(x)^2

% m = (norm(A*(x-c), p)^p-r1^p)/(r2^p-r1^p);
% m = (r1 - tanh(5*x(2))*x(1) - 0.5*(x(2)^2)/(a))/(r1-r2);
% m = (r1 - sign(x(2))*x(1) - 0.5*(x(2)^2)/(a))/(r1-r2);
% m = (r1 - (2*formStep(0.5,1.5,-1,x(2),1)-1)*x(1) - 0.5*(x(2)^2)/(a))/(r1-r2);



m = (r1 - x(1) - sign(x(2))*0.5*(x(2)^2)/(a))/(r1-r2);

f = piecewise(m<=0, 0, exp(-1/m));
f_shift = piecewise(1-m <= 0, 0, exp(-1/(1-m)));

% Normalize the symmetric step function
g = f/(f+f_shift);
bump = g; % Emphasis that g is bump function

if vpa_enable
    bump = vpa(bump,5);
end


end