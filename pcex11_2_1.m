% pcex11_2_1

% Generate simulations of the Ito process X and its approximations using
% explicit order 1.5 strong scheme (2.1)-(2.3)
% 27.02.2015
% =========================================================================
close all
global IFUNC
IFUNC = 2;
nintervals = [8    16    32    64];
Delta = 1./nintervals;
X0 = 1.0;
M = 20; N = 100;
n = length(Delta);
CITaylor = zeros(2,n);
figure
s = idumGenerator(-10,n,100);
errTaylor = zeros(M,1);
for j = 1:n       
    X = zeros(nintervals(j)+1,N);    % matrix to store the exact solution
    Yt = zeros(nintervals(j)+1,N);   % matrix to store the Taylor approximations
    seedsM = idumGenerator(s(j),M,100);
    t = linspace(0,1,nintervals(j)+1);
    for i = 1:M
        seedsN = idumGenerator(seedsM(i),N,100);
        for k = 1:N
            u = gauss_box(seedsN(k),nintervals(j)*2);
            dW = u(1:nintervals(j)) * sqrt(Delta(j));
            dZ = 1/2 * Delta(j)^(3/2) * (u(1:nintervals(j)) + 1/sqrt(3) * u(nintervals(j)+1:end));
            X(:,k) = exactItoSoln(X0,t,dW);
            Yt(:,k) = ExplicitPlaten11_2_1(X0,0,1,dW,dZ);
        end        
        errTaylor(i) = mean(abs(X(end,:) - Yt(end,:)));
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
s = '90% confidence interval for mean approximation error using Platen 11.2.1';
sw = textwrap({s},60);
title(sw);
set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
saveas(h,'pcex11_2_1_ci.jpeg')

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean mean approximation error using Platen 11.2.1 against time step size \Delta in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex11_2_1_log.jpeg')

