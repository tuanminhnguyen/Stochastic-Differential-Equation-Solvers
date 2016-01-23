% Compute the Euler-Maruyama approximations to the solution of an SDE
% See Kloden, chapter 9.

%            X0     initial condition
%            a      drift constant
%            b      diffusion constant
% Inputs     Delta  constant time step size
%            dW     vector storing Wiener increments

% Output     Y      vector storing exact solution evaluated at each grid
%                   grid point in t
% 27.02.2015
%==========================================================================
function Y = EulerApproxSDE(X0,t0,tf,dW)
global IFUNC
if IFUNC == 1
    a = 1.5;
    b = 1.0;
else
    a = 1.5;
    b = 0.1;
end
m = length(dW);
Delta = (tf - t0) / m;
aDelta = a * Delta;
Y = zeros(m+1,1);
Y(1) = X0;

for n = 1:m
    Y(n+1) = (1 + aDelta + b*dW(n)) * Y(n); % chapter 9 (2.5)
end

end