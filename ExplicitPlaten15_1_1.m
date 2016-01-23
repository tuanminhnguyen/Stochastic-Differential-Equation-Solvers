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
function Y = ExplicitPlaten15_1_1(X0,t0,tf,dW)
m = length(dW);
Delta = (tf - t0) / m;
Y = zeros(m+1,1);
Y(1) = X0;
dWsq = dW.^2;
t = linspace(t0,tf,m+1);
for n = 1:m
    [a,b,ax,bx] = func(Y(n),t(n));
    [gamp,gamm,gam] = func_gam(Y(n),ax,bx,Delta,dW(n));
    [~,~,agam,~] = func(gam,t(n));
    [~,~,~,bgamp] = func(gamp,t(n));
    [~,~,~,bgamm] = func(gamm,t(n));
    Y(n+1) = Y(n) + 1/2 * (agam + ax) * Delta + ...
        1/4 * (bgamp + bgamm + 2*bx) * dW(n) + ...
        1/4 * (bgamp - bgamm) * (dWsq(n) - Delta) * Delta^(-1/2);
end
end

function [a,b,ax,bx] = func(x,t)
a = 1.5;
b = 0.01;
ax = a * x;
bx = b * x;
end

function [gamp,gamm,gam] = func_gam(y,a,b,Delta,dW)
gamp = y + a * Delta + b * sqrt(Delta);
gamm = y + a * Delta - b * sqrt(Delta);
gam = y + a * Delta + b * dW;
end