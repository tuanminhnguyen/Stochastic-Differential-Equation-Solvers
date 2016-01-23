close all
global IFUNC
IFUNC = 3;
nnodes = [32 64 128 256];  %[128 256 512 1024];
M = 20; N = 100;
X0 = 0.1;
CI = zeros(2,4);
err = zeros(M,1);
for j = 1:4    
    X = zeros(nnodes(j) + 1,N);
    Y = zeros(nnodes(j) + 1,N);
    t = linspace(0,1,nnodes(j)+1);
    for i = 1:M
        for k = 1:N
            dW = sqrt(1/nnodes(j)).*randn(nnodes(j),1) + 0;        
            X(:,k) = exactItoSoln(X0,t,dW);
            Y(:,k) = ExplicitPlaten15_1_1(X0,0,1,dW);        
        end              
        err(i) = mean(Y(end,:)) - mean(X(end,:));
    end
     [~,~,CI(:,j),~] = normfit(err,0.1);   
end

Delta = 1./nnodes;

halflength = (CI(2,:) - CI(1,:))/2;
midpoint = (CI(2,:) + CI(1,:))/2;
% h = errorbar(Delta, midpoint, halflength,'LineStyle','none');
% xlabel('$\Delta$ = time step size','interpreter','LaTex')
% ylabel('$\epsilon$ = average approximation error','interpreter','LaTex')
% s = '90% confidence interval for mean approximation error using Weak Euler Approximation';
% sw = textwrap({s},60);
% title(sw);
% set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
% saveas(h,'9_3_3_CIapproxErr.jpeg')

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean error \epsilon against time step size |\Delta| in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex15_1_1.jpeg')
