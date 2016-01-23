% Compute the explitcit order 1.0 strong approximation to the solution of an SDE
% See Kloden, 11.1 (1.3) pg 374.

% Inputs     X0     initial condition
%            t0     starting time
%            tf     final time
%            dW     vector storing Wiener increments
%            dZ     vector storing Wiener increments

% Output     Y      vector storing exact solution evaluated at each grid
%                   grid point in t
% 27.02.2015
%==========================================================================
function Y = ExplicitPlaten11_2_1(X0,t0,tf,dW,dZ)
m = length(dW);
Delta = (tf - t0) / m;
Y = zeros(m+1,1);
Y(1) = X0;
dWsq = dW.^2;
t = linspace(t0,tf,m+1);
for n = 1:m
    [a,b,ax,bx] = func(Y(n),t(n));
    [gamp,gamm] = func_gam(Y(n),ax,bx,Delta);
    [~,~,agamp,bgamp] = func(gamp,t(n));
    [~,~,agamm,bgamm] = func(gamm,t(n));
    [phip,phim] = func_phi(gamp,bgamp,Delta);
    [~,~,~,bphip] = func(phip,t(n));
    [~,~,~,bphim] = func(phim,t(n));
    Y(n+1) = Y(n) + bx * dW(n) + 1/2/sqrt(Delta) * (agamp - agamm) * dZ(n) +...
        1/4 * (agamp + 2*ax + agamm) * Delta + ...
        1/4/sqrt(Delta) * (bgamp - bgamm) * (dWsq(n) - Delta) + ...
        1/2/Delta * (bgamp - 2*bx + bgamm) * (dW(n) * Delta - dZ(n)) + ...
        1/4/Delta * (bphip - bphim - bgamp + bgamm) * (1/3*dWsq(n) - Delta) * dW(n);    
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

function [gamp,gamm] = func_gam(y,a,b,Delta)
gamp = y + a * Delta + b * sqrt(Delta);
gamm = y + a * Delta - b * sqrt(Delta);
end

function [phip,phim] = func_phi(gamp,bgamp,Delta)
phip = gamp + bgamp * sqrt(Delta);
phim = gamp - bgamp * sqrt(Delta);
end