% Compute the Euler-Maruyama approximations to the solution of an SDE
% See Kloden, 10.2 (2.1) pg 340.

% Inputs     X0     initial condition
%            t0     starting time
%            tf     final time
%            dW     vector storing Wiener increments

% Output     Y      vector storing exact solution evaluated at each grid
%                   grid point in t
% 27.02.2015
%==========================================================================
function Y = MilsteinApproxSDE(X0,t0,tf,dW)
m = length(dW);
Delta = (tf - t0) / m;
Y = zeros(m+1,1);
Y(1) = X0;
dWsq = dW.^2;
t = linspace(t0,tf,m);
for n = 1:m 
    [a,b,bp] = func(Y(n),t(n));
    Y(n+1) = Y(n) + a * Delta + b * dW(n) + b * bp * (dWsq(n) - Delta) / 2; % chapter 9 (2.5)
end
end

function [a,b,bp] = func(x,t)
global IFUNC
if IFUNC == 1    
    a = 1.5;
    b = 1.0;
else
    a = 1.5;
    b = 0.1;
end
a = a * x;
b = b * x;
bp = b;
end