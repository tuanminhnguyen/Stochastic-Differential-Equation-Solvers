nnodes = [4 8 16 32];  %[128 256 512 1024];
exact = 0.1 * exp(3/2);

M = 20; N = 100;
X0 = 0.1;
CI = zeros(2,4);
err = zeros(M,1);
for j = 1:4
    Y = zeros(nnodes(j) + 1,N);
    Delta = 1/nnodes(j);
    sqrtDelta = sqrt(1/nnodes(j));    
    for i = 1:M
        for k = 1:N
            dW1 = sqrtDelta.*randn(nnodes(j),1) + 0;
            dW2 = sqrtDelta.*randn(nnodes(j),1) + 0;
            Ym = zeros(nnodes(j) + 1,1);
            Ym(1,1) = 0.1;
            for n = 1:nnodes(j)
                Ym(n+1,1) = (1 + 1.5 * Delta + 0.1*(dW1(n) + dW2(n)))*Ym(n);
            end
            Y(:,k) = Ym;
        end
        
        err(i) = mean((Y(end,:))) - exact;
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
s = 'plot of mean Euler error \epsilon against time step size |\Delta| in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex14_1_2.jpeg')
