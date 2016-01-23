% pcex9_4_1

close all
global IFUNC
IFUNC = 2;
a = 1.5; b = 0.1; X0 = 1.0;
% compute mean simulation errors
m = 100;
ninterval = 16;
N = 100;
% Compute the mean of the Ito process:
%     E(X_T) = E(X_0) exp(aT)
meanIto = X0 * exp(a * 1);
Y = zeros(ninterval+1,N);            % matrix to store the approximate sol   
s = idumGenerator(-2,N,100);
hat_mu = zeros(m,1); 

for j = 1:m
    seedsN = idumGenerator(s(j),N,100);   
    for i = 1:N
        dW = WienerIncrement(0,1,ninterval,seedsN(i));
        Y(:,i) = eulerMaruyama(X0,ninterval,dW);           
    end       
    hat_mu(j) = sum(Y(end,:) - meanIto ) / N; % chapter 9 (4.4)
end

% compute confidence intervals
M = [10 20 40 100];
CI = zeros(2,4);

for j = 1:4    
    [~,~,CI(:,j),~] = normfit(hat_mu(1:M(j)),0.1);
end
halflength = (CI(2,:) - CI(1,:))/2;
midpoint = (CI(2,:) + CI(1,:))/2;
h = errorbar(M, midpoint, halflength,'LineStyle','none');
xlabel('M = number of simulation batches')
ylabel('$\hat \mu$ = mean approximation error','interpreter','LaTex')
r = '90% confidence interval for mean error estimate with varying number of batches';
sw = textwrap({r},60);
title(sw);
set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
saveas(h,'pcex9_4_1.jpeg')