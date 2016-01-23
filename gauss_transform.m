% Transform a pair of independent standard Gaussian random values (x1,x2)
% into a pair of random values drawn from two independent Gaussian
% distributions with covariance matrix C and mean vector mu
% ( E[X1] = mu_1,  E[X2] = mu_2 ).

% This function uses the transform method in Kloeden & Patten, 19.

% Inputs  x   2xn matrix storing n pairs of ind. standard Gaussian random values
%         mu  2x1 array storing the means of the two new random variables 
%         C   covariance matrix for the new random variables
%         n   optional number of pairs desired, DEFAULT = 1

% Output  x   2xn

% 21.02.2015
% =========================================================================
function x = gauss_transform(x,mu,C)
% inverting the covariance matrix and form lower Cholesky matrix
L = chol(C,'lower')'; % the transpose is needed for the left multiplication below
x =  x * L;
x(:,1) = x(:,1) + mu(1,1);
x(:,2) = x(:,2) + mu(2,1);
end
