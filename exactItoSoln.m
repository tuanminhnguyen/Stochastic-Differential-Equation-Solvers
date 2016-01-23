% Computes the exact solution to the SDE
%    dX_t = aX_t dt + bX_t dW_t
% See chapter 9 (2.5)

%            X0     initial condition
% Inputs     t      vector storing the partition grid points
%            dW     vector storing Wiener increments

% Output     X      vector storing exact solution evaluated at each grid
%                   grid point in t
%==========================================================================
function X = exactItoSoln(X0,t,dW)
global IFUNC
if IFUNC == 1
    a = 1.5;
    b = 1.0;
    X0 = 1.0;
elseif IFUNC == 2
    a = 1.5;
    b = 0.1;
else
    a = 1.5;
    b= 0.01;
end

sumdW = 0;
m = length(dW);
X = zeros(m+1,1);
X(1) = X0;

for n = 1:m
    sumdW = sumdW + dW(n); % dW_0 = 0 by definition of Wiener process
    X(n+1) = exp((a-1/2*b^2) * t(n+1) + b * sumdW);
end
X = X0 * X;
end
