% pcex10_4_2

% compute the 90% Confidence Interval for the mean error approximations using
% 1.5 order strong Taylor approximation, for different Delta's
% 27.02.2015
% =========================================================================
close all
global IFUNC
IFUNC = 2;
nintervals = [8    16    32    64];
X0 = 1.0;
M = 20; N = 100;
n = length(nintervals);
CITaylor = zeros(2,n);
figure
s = idumGenerator(-10,n,100);
Delta = 1./nintervals;
for j = 1:n       
    X = zeros(nintervals(j)+1,N);    % matrix to store the exact solution
    Yt = zeros(nintervals(j)+1,N);   % matrix to store the Taylor approximations
    seedsM = idumGenerator(s(j),M,100);
    t = linspace(0,1,nintervals(j)+1);
    errTaylor = zeros(M,1);    
    for i = 1:M
        seedsN = idumGenerator(seedsM(i),N,100);
        for k = 1:N
            u = gauss_box(seedsN(k),nintervals(j)*2);
            dW = u(1:nintervals(j)) * sqrt(Delta(j));
            dZ = 1/2 * Delta(j)^(3/2) * (u(1:nintervals(j)) + 1/sqrt(3) * u(nintervals(j)+1:end));
            X(:,j) = exactItoSoln(X0,t,dW);
            Yt(:,j) = TaylorApproxSDE(X0,0,1,dW,dZ);
        end        
        errTaylor(i) = mean(sum(abs(X(end,:) - Yt(end,:))));
    end
    [~,~,CITaylor(:,j),~] = normfit(errTaylor,0.1);
end

%%
% plot Euler

halflength = (CITaylor(2,:) - CITaylor(1,:))/2;
midpoint = (CITaylor(2,:) + CITaylor(1,:))/2;
h = errorbar(Delta, midpoint, halflength,'LineStyle','none');
xlabel('$\Delta$ = time step size','interpreter','LaTex')
ylabel('$\epsilon$ = average approximation error','interpreter','LaTex')
s = '90% confidence interval for Taylor mean approximation error';
sw = textwrap({s},60);
title(sw);
set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
saveas(h,'pcex10_4_2_ci.jpeg')

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean Taylor error \epsilon against time step size \Delta in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex10_4_2_log.jpeg')

