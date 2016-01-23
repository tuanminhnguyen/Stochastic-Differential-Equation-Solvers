% Compute the Taylor approximations to the solution of an SDE
% See Kloden, chapter 10.

% Inputs     X0     initial condition
%            t0     starting time
%            tf     final time
%            dW     Wiener increments
%            dZ     Wiener increments

% Output     Y      vector storing exact solution evaluated at each grid
%                   grid point in t
% 27.02.2015
%==========================================================================
function Y = TaylorApproxSDE(X0,t0,tf,dW,dZ)
m = length(dW);
Delta = (tf - t0) / m;
t = linspace(t0,tf,m+1);
Y = zeros(m+1,1);
Y(1) = X0;
dWsq = dW.^2;
for n = 1:m
    [a,b,ap,bp,app,bpp] = func(Y(n),t(n));
    Y(n+1) = Y(n) + a*Delta + b*dW(n) + 1/2 * b * bp * (dWsq(n) - Delta) + ...
        ap * b * dZ(n) + 1/2 *(a*ap + 1/2 * b^2 * app) * Delta^2 + ...
        (a * bp + 1/2 * b^2 * bpp) * (dW(n) * Delta - dZ(n)) + ...
        1/2 * b * (b * bpp + bp^2) * ( 1/3 * dWsq(n) - Delta) * dW(n); % chapter 10 (4.1 pg 351)
end
end


function [a,b,ap,bp,app,bpp] = func(x,t)
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
ap = a;
bp = b;
app = 0;
bpp = 0;
end