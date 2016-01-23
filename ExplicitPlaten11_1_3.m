% Compute the explitcit order 1.0 strong approximation to the solution of an SDE
% See Kloden, 11.1 (1.3) pg 374.

% Inputs     X0     initial condition
%            t0     starting time
%            tf     final time
%            dW     vector storing Wiener increments

% Output     Y      vector storing exact solution evaluated at each grid
%                   grid point in t
% 27.02.2015
%==========================================================================
function Y = ExplicitPlaten11_1_3(X0,t0,tf,dW)
m = length(dW);
Delta = (tf - t0) / m;
Y = zeros(m+1,1);
Y(1) = X0;
dWsq = dW.^2;
t = linspace(t0,tf,m+1);
for n = 1:m
    [a,b,ax,bx] = func(Y(n),t(n));
    ybar = func_y(Y(n),a,b,Delta);
    [~,~,~,bbar] = func(ybar,t(n));
    Y(n+1) = Y(n) + ax * Delta + bx * dW(n) + 1/2/sqrt(Delta) *  ...
        (bbar - bx) * (dWsq(n) - Delta); % chapter 11 (1.3)
end
end

function [a,b,ax,bx] = func(x,t)
global IFUNC
if IFUNC == 1
    a = 1.5;
    b = 1.0;
else
    a = 1.5;
    b = 0.1;
end
ax = a * x;
bx = b * x;
end

function ybar = func_y(y,a,b,Delta)
ybar = y + a * Delta + b * sqrt(Delta);
end