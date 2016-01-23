% pcex10_3_2

% compute the 90% Confidence Interval for the mean error approximations using
% Euler-Maruyama and Milstein schemes, for different Delta's
% 27.02.2015
% =========================================================================
close all
global IFUNC
IFUNC = 2;
nintervals = [8    16    32    64];
X0 = 1.0;
M = 20; N = 100;
n = length(nintervals);
CIEuler = zeros(2,4);
CIMilstein = zeros(2,4);
s = idumGenerator(-2,n,100);
for j = 1:n
    X = zeros(nintervals(j)+1,N);    % matrix to store the exact solution
    Ye = zeros(nintervals(j)+1,N);   % matrix to store the Milstein approximations
    Ym = zeros(nintervals(j)+1,N);   % matrix to store the Euler approximations
    errEuler = zeros(M,1);    
    errMilstein = zeros(M,1);    
    seedsM = idumGenerator(s(j),M,100);
    t = linspace(0,1,nintervals(j)+1);
    for i = 1:M
        seedsN = idumGenerator(seedsM(i),N,100);
        for k = 1:N
            dW = WienerIncrement(0,1,nintervals(j),seedsN(k));
            X(:,j) = exactItoSoln(X0,t,dW);
            Ye(:,j) = EulerApproxSDE(X0,0,1,dW);
            Ym(:,j) = MilsteinApproxSDE(X0,0,1,dW);
        end        
        errEuler(i) = mean(sum(abs(X(end,:) - Ye(end,:))));
        errMilstein(i) = mean(sum(abs(X(end,:) - Ym(end,:))));
    end
    [~,~,CIEuler(:,j),~] = normfit(errEuler,0.1);
    [~,~,CIMilstein(:,j),~] = normfit(errMilstein,0.1);    
end


%%
% % plot Euler
Delta = 1./nintervals;
halflength = (CIEuler(2,:) - CIEuler(1,:))/2;
midpoint = (CIEuler(2,:) + CIEuler(1,:))/2;
% h = errorbar(Delta, midpoint, halflength,'LineStyle','none');
% xlabel('$\Delta$ = time step size','interpreter','LaTex')
% ylabel('$\epsilon$ = average approximation error','interpreter','LaTex')
% s = '90% confidence interval for Euler mean approximation error';
% sw = textwrap({s},60);
% title(sw);
% set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
% saveas(h,'pcex10_3_2_euler.jpeg')

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean Euler error \epsilon against time step size \Delta in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex10_3_2_euler.jpeg')


% figure
% plot Milstein
% halflength = (CIMilstein(2,:) - CIMilstein(1,:))/2;
% midpoint = (CIMilstein(2,:) + CIMilstein(1,:))/2;
% h = errorbar(Delta, midpoint, halflength,'LineStyle','none');
% xlabel('$\Delta$ = time step size','interpreter','LaTex')
% ylabel('$\epsilon$ = average approximation error','interpreter','LaTex')
% s = '90% confidence interval for Milstein mean approximation error';
% sw = textwrap({s},60);
% title(sw);
% set(gca, 'FontName', 'Latin Modern Roman', 'FontSize', 13);
% saveas(h,'9_3_3_CIapproxErr.jpeg')

figure
h = plot(log2(Delta), log2(abs(midpoint)));
xlabel('$$\log_2 \Delta$$','interpreter','LaTex')
ylabel('$$\log_2 \epsilon$$','interpreter','LaTex')
s = 'plot of mean Milstein error \epsilon against time step size \Delta in log2 scale';
sw = textwrap({s},60);
title(sw,'interpreter','LaTex');
set(gca, 'FontSize', 13)
saveas(h,'pcex10_3_2_milstein.jpeg')
