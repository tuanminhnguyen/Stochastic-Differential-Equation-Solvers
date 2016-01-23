% pcex9_4_2, 9_4_3
close all
M = 20;
N = 100;
nintervals = [8    16    32    64];
% delta = (2 * ones(4,1)) .^ ((-3:-1:-6)');
n = length(nintervals);

global IFUNC
IFUNC = 2;
X0 = 1.0; a = 1.5; b = 0.1;
% Compute the mean of the Ito process
%     E(X_T) = E(X_0) exp(aT)
meanIto = X0 * exp(a * 1);
s = idumGenerator(-2,n,100);
hat_mu = zeros(M,1);
mu = zeros(4,1);
CI = zeros(2,4);
for j = 1:n    
    seedsK = idumGenerator(s(j),M,100);
    for k = 1:M       
        Y = zeros(nintervals(j)+1,N);            % matrix to store the approximate soln
        seedsN = idumGenerator(seedsK(k),N,100);
        for i = 1:N
            dW = WienerIncrement(0,1,nintervals(j),seedsN(i));
            Y(:,i) = eulerMaruyama(X0,nintervals(j),dW);
        end
        hat_mu(k) = abs(sum(Y(end,:) - meanIto )) / N; % chapter 9 (4.4)
    end
    mu(j) = mean(hat_mu); % chapter 9 (4.5)
    [~,~,CI(:,j),~] = normfit(hat_mu,0.1);
end

% compute confidence intervals

halflength = (CI(2,:) - CI(1,:))/2;
midpoint = (CI(2,:) + CI(1,:))/2;
delta = 1./nintervals;
h1 = errorbar(delta, midpoint, halflength,'LineStyle','none');
xlabel('$$\Delta$$ = time step size','interpreter','LaTex');
ylabel('$$\hat \mu$$ = mean approximation error','interpreter','LaTex')
s = '90% confidence interval for mean error estimate with varying time step size';
sw = textwrap({s},60);
set(gca, 'FontSize', 13)
title(sw);

saveas(h1,'pcex9_4_3ci.jpeg')

% pcex9_4_3
figure
h = plot(log2(delta), log2(abs(mu)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \hat \mu$$','interpreter','LaTex')
s = 'plot of mean_mu against time step size \Delta in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex9_4_3log.jpeg')