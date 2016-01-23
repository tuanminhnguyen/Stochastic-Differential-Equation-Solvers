% pcex9.3.1
% Estimate the Variance of the absolute error measure e (chapter 9 (3.1))
% for Euler pathwise approximation, use this variance to construct a
% Confidence Interval for e.
% See chapter 9 (3.3), (3.4), (3.5) pg 312.

% 26.02.2015
%==========================================================================
close all
global IFUNC
IFUNC = 2;
% compute simulation errors
m = 100; % number of simulation batches
err = zeros(m,1);
s = idumGenerator(-50,m,100);
ninterval = 16;
X0 = 1.0;
t = linspace(0,1,ninterval+1);
N = 100;        % number of trajectories simulated for each batch
X = zeros(ninterval+1,N);            % matrix to store the exact solution
Y = zeros(ninterval+1,N);
for i = 1:m
    seedsN = idumGenerator(s(i),N,100);
    for j = 1:N
        dW = WienerIncrement(0,1,ninterval,seedsN(j));
        X(:,j) = exactItoSoln(X0,t,dW);
        Y(:,j) = eulerMaruyama(X0,ninterval,dW);
    end    
    Sum = sum(abs(X(end,:) - Y(end,:)));
    err(i) = Sum/N;
end

% compute confidence intervals
M = [10 20 40 100];
CI = zeros(2,4);

for j = 1:4
    [~,~,CI(:,j),~] = normfit(err(1:M(j)),0.1);
end

% plot confidence intervals
figure
hold on
halflength = (CI(2,:) - CI(1,:))/2;
midpoint = (CI(2,:) + CI(1,:))/2;
h = errorbar(M, midpoint, halflength,'LineStyle','none');
xlabel('M = number of simulation batches','interpreter','LaTex')
ylabel('$\epsilon$ = average approximation error','interpreter','LaTex')
r = '90% confidence interval for mean approximation error with varying number of batches';
sw = textwrap({r},60);
title(sw);
set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
saveas(h,'pcex9_3_3.jpeg')
