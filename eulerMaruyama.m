% Compute the Euler-Maruyama approximations to the solution of an SDE
% See Kloden, chapter 9.

% Inputs     X0     initial condition
%            nintervals  number of partition intervals on domain
%            dW     vector storing Wiener increments

% Output     Y      vector storing exact solution evaluated at each grid
%                   grid point in t
% 27.02.2015
%==========================================================================
function Y = eulerMaruyama(X0,nintervals,dW)
global IFUNC
if IFUNC == 1
    a = 1.5;
    b = 1.0;
elseif IFUNC == 2
    a = 1.5;
    b = 0.1;
else
    a = 1.5;
    b= 0.01;
end
aDelta = a / nintervals; % = a * Delta
m = length(dW);
Y = zeros(m+1,1);
Y(1) = X0;

for n = 1:m
    Y(n+1) = (1 + aDelta + b*dW(n)) * Y(n); % chapter 9 (2.5)
end

end