X0 = [1; 0];
global a b A B
a = 5; b = 0.08; T = 1;
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
             dW = randn(ninterv(j),1);
             X = exactItoSoln2D(X0,t,dW)';
             
             for n = 1:ninterv(j)
                gamp = (eye(2) + A * Delta + B * sqrt(Delta)) * Y(:,n);
                gamm = (eye(2) + A * Delta - B * sqrt(Delta)) * Y(:,n);
                gam = (eye(2) +  A * Delta + B * dW(n)) * Y(:,n);
                phi = 1/4 * B * (gamp + gamm + 2 * Y(:,n)) * dW(n) + ...
                      1/4 * B * (gamp - gamm) * (dW(n)^2 - Delta) * Delta^(-1/2);
                Ybar = Y(:,n) + 1/2 * A * (gam + Y(:,n)) * Delta + phi;
                Y(:,n+1) = Y(:,n) + 1/2 * A * (Ybar + Y(:,n)) * Delta + phi;
             end
             Yk(:,k) = Y(:,end);
             Xk(:,k) = X(:,end);
        end
        YN = Yk(1,:).^2 + Yk(2,:).^2; % I'm not sure how to measure the error here
        % anyway the result is clearly incorrect.
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
saveas(h,'pcex15_5_3.jpeg')
