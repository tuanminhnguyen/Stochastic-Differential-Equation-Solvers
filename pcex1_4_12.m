% Random number generator for a pair of normal random numbers with
% mean vector mu and (co)variances E(X1^2) = h, E(X2^2) = h^3/3, E(X1X2) = h^2
% for any h > 0.

% Inputs    n    number of pairs of values
%           h    desired parameter for variance and covariance
%           mu   2x1 mean vector

% Output    x    2xn matrix storing pairs of Gaussian random values satisfy
%                the above conditions

% 21.02.2015
% =========================================================================
function x = pcex1_4_12(n,mu,h)
% form covariance matrix C
C = [ h     ,  h^2/2;...
      h^2/2 ,  h^3/3];
  
% produce pair(s) of standard Gaussian random numbers using Box-Muller or Polar
% Marsaglia
if nargin == 2
    n = 1;
end
x = gauss_box(-2, n);
x = reshape(x,n/2,2);

% transform into the new Gaussian pairs
x = gauss_transform(x,mu,C);
end