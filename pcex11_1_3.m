% pcex11_1_3

% Generate simulations of the Ito process X and its approximations
% 27.02.2015
% =========================================================================
close all
global IFUNC
IFUNC = 2;
nintervals = [8    16    32    64];
Delta = 1./nintervals;
X0 = 1.0;
M = 20; N = 100;
n = length(nintervals);
CITaylor = zeros(2,n);
s = idumGenerator(-10,n,100);
errPlaten = zeros(M,1);
for j = 1:n       
    X = zeros(nintervals(j)+1,N);    % matrix to store the exact solution
    Yt = zeros(nintervals(j)+1,N);   % matrix to store the Taylor approximations
    seedsM = idumGenerator(s(j),M,100);
    t = linspace(0,1,nintervals(j)+1);
    for i = 1:M
        seedsN = idumGenerator(seedsM(i),N,100);
        for k = 1:N            
            dW = WienerIncrement(0,1,nintervals(j),seedsN(k));
            X(:,j) = exactItoSoln(X0,t,dW);
            Yt(:,j) = ExplicitPlaten11_1_3(X0,0,1,dW);
        end        
        errPlaten(i) = mean(sum(abs(X(end,:) - Yt(end,:))));
    end
    [~,~,CITaylor(:,j),~] = normfit(errPlaten,0.1);
end

%%
% plot Euler

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean approximation error against time step size \Delta in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex11_1_3.jpeg')

