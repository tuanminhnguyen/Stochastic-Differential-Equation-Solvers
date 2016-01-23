% pcex9_2_1
% plot a sample of equidistant Euler approximations and corresponding
% exact solutions
global IFUNC
IFUNC = 1;
n = 10;
nintervals = 32;
figure
delta = 1/nintervals;
X0 = 1.0;
t = linspace(0,1,nintervals+1);

dW = reshape(WienerIncrement(0,1,nintervals*n,-2),n,nintervals);
for i = 1:n
    % Euler approximations
    Y = eulerMaruyama(X0,nintervals,dW(i,:));
    % exact solution
    X = exactItoSoln(X0,t,dW);
    plot(t, Y, 'ro-', t, X, 'g*-')
    hold on
end

title(sprintf('%d samples of equidistant Euler approximations and exact solutions',n))
legend('approximate solution','exact solution');
set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);

% saveas(h,'pcex9_2_1.jpeg')
hold off