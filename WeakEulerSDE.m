
function Y = WeakEulerSDE(X0,t0,tf,dW1,dW2)
m = length(dW1);
Delta = (tf - t0) / (m - 1);
Y = zeros(m,1);
Y(1) = X0;
t = linspace(t0,tf,m);
for n = 1:m-1
    [~,~,~,ax,bx1,bx2] = func(Y(n),t(n));
    Y(n+1) = Y(n) + ax * Delta + bx1 * dW1(n) + bx2 * dW2(n);
%     Y(n+1,1) = (1 + 1.5 * Delta + 0.1*(dW1(n) + dW2(n)))*Y(n);
end
% Y(end-3:end)
end

function [a,b1,b2,ax,bx1,bx2] = func(x,t)
a = 1.5;
b1 = 0.1;
b2 = 0.1;
ax = a * x;
bx1 = b1 * x;
bx2 = b2 * x;
end