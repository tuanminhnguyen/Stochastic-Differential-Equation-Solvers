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
function Y = ImplicitEuler(X0,t0,tf,dW)






m = length(dW);
Delta = (tf - t0) / (m - 1);
Y = zeros(m,1);
Y(1) = X0;
dWsq = dW.^2;
t = linspace(t0,tf,m);
for n = 1:m-1
    [a,b,ax,bx] = func(Y(n),t(n));
    [gamp,gamm] = func_gam(Y(n),a,b,Delta);
    [~,~,agamp,bgamp] = func(gamp,t(n));
    [~,~,agamm,bgamm] = func(gamm,t(n));
    [phip,phim] = func_phi(gamp,bgamp,Delta);
    [~,~,~,bphip] = func(phip,t(n));
    [~,~,~,bphim] = func(phim,t(n));
    Y(1,n+1) = Y(1,n) + (alpha1 * dW(n) + 1/2/sqrt(Delta) * (agamp - agamm) * dZ(n) +...
        1/4 * (agamp + 2*ax + agamm) * Delta + ...
        1/4/sqrt(Delta) * (bgamp - bgamm) * (dWsq(n) - Delta) + ...
        1/2/Delta * (bgamp - 2*bx + bgamm) * (dW(n) * Delta - dZ(n)) + ...
        1/4/Delta * (bphip - bphim - bgamp + bgamm) * (1/3*dWsq(n) - Delta) * dW(n);
    
end
end
