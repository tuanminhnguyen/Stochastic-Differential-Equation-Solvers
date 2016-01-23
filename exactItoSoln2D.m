% Computes the exact solution to the SDE
%    dX_t = aX_t dt + bX_t dW_t
% See chapter 9 (2.5)

%            X0     initial condition
% Inputs     t      vector storing the partition grid points
%            dW     vector storing Wiener increments

% Output     X      vector storing exact solution evaluated at each 
%                   grid point in t
%==========================================================================
function X = exactItoSoln2D(X0,t,dW)
global A B
a = 5; b = 0.01;
% A = [-a a; a -a]; B = [b 0; 0 b];
sumdW = 0;
m = length(dW);
X = zeros(m+1,2);
X(1,:) = X0;
P = 1/sqrt(2) *[1 1;1 -1];

for n = 1:m
    sumdW = sumdW + dW(n); % dW_0 = 0 by definition of Wiener process
%     X(n+1,:) = exp((A - 1/2 * B^2) * t(n+1) + B * sumdW) * X0;    
    X(n+1,:) = P*[exp((-a-1/2*b^2+a)*t(n+1) + b*sumdW), 0;...
                    0, exp((-a-1/2*b^2-a)*t(n+1) + b*sumdW)] * P * X0;
end
% X = X * X0;
end
