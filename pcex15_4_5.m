X0 = [1; 0];
global a b A B
a = 5; b = 0.01; T = 1;
A = [-a a;a -a]; B = [b 0; 0 b];Bs = [b^2 0 ; 0 b^2];
M = 20; N = 100;
m = 4; ninterv =  (2 * ones(1,m)).^(2:5);

CI = zeros(2,m);
err = zeros(M,1);

for j = 1:m
    nnodes = ninterv(j) + 1;    
    Delta = 1/ninterv(j);
    Yk = zeros(2,N);
    Xk = zeros(2,N);
    Y = zeros(2,nnodes);    
    Y(:,1) = X0;
    t = linspace(0,T,nnodes);
    for i = 1:M
        for k = 1:N
             dW = rand(ninterv(j),1);
             dW(dW>0.5) = sqrt(Delta);
             dW(dW<=0.5) = -sqrt(Delta);
             X = exactItoSoln2D(X0,t,dW)';
             
             for n = 1:ninterv(j)
                Y(:,n+1) = (eye(2) - (A - Bs*Delta) - B * dW(n)) \ Y(:,n);
             end
             Yk(:,k) = Y(:,end);
             Xk(:,k) = X(:,end);
        end
        YN = Yk(1,:).^2 + Yk(2,:).^2;
        XN = Xk(1,:).^2 + Xk(2,:).^2;
        err(i) = mean(YN) - mean(XN);
    end
    [~,~,CI(:,j),~] = normfit(err,0.1);
end
             
Delta = 1./ninterv;

halflength = (CI(2,:) - CI(1,:))/2;
midpoint = (CI(2,:) + CI(1,:))/2;
% h = errorbar(Delta, midpoint, halflength,'LineStyle','none');
% xlabel('$\Delta$ = time step size','interpreter','LaTex')
% ylabel('$\epsilon$ = average approximation error','interpreter','LaTex')
% s = '90% confidence interval for mean approximation error';
% sw = textwrap({s},60);
% title(sw);
% set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean error \epsilon against time step size |\Delta| in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex15_4_5.jpeg')
