% test pcex1_4_12
% 21.02.2015
% =========================================================================

% experiment 1
n = 1E4;
mu= [0;0];
h = 1;
x = pcex1_4_12(n,mu,h);

figure
subplot(2,1,1)
[mu, var] = chk_gauss_normal(x(:,1)) %#ok<*ASGLU,*NOPTS>
subplot(2,1,2)
[mu, var] = chk_gauss_normal(x(:,2))

cov(x(:,1),x(:,2))


clear

% experiment 2
n = 1E4;
mu= [1;2];
h = 2;
x = pcex1_4_12(n,mu,h);

figure
subplot(2,1,1)
[mu, var] = chk_gauss_normal(x(:,1)) %#ok<*ASGLU,*NOPTS>
subplot(2,1,2)
[mu, var] = chk_gauss_normal(x(:,2))

cov(x(:,1),x(:,2))